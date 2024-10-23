#! /usr/bin/env python
import pysam  # Import the pysam module
import argparse  # Import the argparse module
from collections import defaultdict  # Import the defaultdict class from the collections module
import sys  # Import the sys module
import multiprocessing


def get_chromosome_ranges(bam_file, num_chunks):
    """
    Divide the BAM file's reference chromosomes into ranges for parallel processing.
    Ensure that each chunk range is exclusive of the end to avoid overlapping processing.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    chromosome_ranges = []

    for chromosome in bam.references:
        length = bam.lengths[bam.references.index(chromosome)]
        chunk_size = length // num_chunks
        for start in range(0, length, chunk_size):
            end = min(start + chunk_size, length)
            chromosome_ranges.append((chromosome, start, end))

    bam.close()
    return chromosome_ranges


def process_bam_chunk(bam_file, variants, min_supporting_variants, continuity_threshold, switch_threshold,
                      dominant_fraction_threshold, min_mapping_quality, min_read_length, chunk, min_variant_quality):
    chromosome, start, end = chunk
    analyzer = HaplotypeAnalyzer(variants, min_supporting_variants, continuity_threshold,
                                 switch_threshold, dominant_fraction_threshold, min_variant_quality, k_consecutive=3)
    results = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chromosome, start, end):
            if (read.is_secondary or read.is_supplementary or read.mapping_quality < min_mapping_quality
                or read.query_length < min_read_length):
                continue
            result = analyzer.calculate_score(read)
            results.append(result)

    return results


def process_bam(bam_file, variants, min_supporting_variants, continuity_threshold, switch_threshold,
                dominant_fraction_threshold, min_mapping_quality, min_read_length, min_variant_quality, num_workers=4):
    chromosome_ranges = get_chromosome_ranges(bam_file, num_workers)

    with multiprocessing.Pool(processes=num_workers) as pool:
        results = pool.starmap(
            process_bam_chunk,
            [(bam_file, variants, min_supporting_variants, continuity_threshold, switch_threshold,
              dominant_fraction_threshold, min_mapping_quality, min_read_length, chunk, min_variant_quality)
             for chunk in chromosome_ranges]
        )

    return [item for sublist in results for item in sublist]



class HaplotypeAnalyzer:
    def __init__(self, variants, min_supporting_variants, continuity_threshold, switch_threshold,
                 dominant_fraction_threshold=0.9, min_variant_quality=20, k_consecutive=3):
        self.variants = variants
        self.min_supporting_variants = min_supporting_variants
        self.continuity_threshold = continuity_threshold
        self.switch_threshold = switch_threshold
        self.dominant_fraction_threshold = dominant_fraction_threshold
        self.previous_read_info = {}  # Dictionary to store the last read's information
        self.min_variant_quality = min_variant_quality
        self.k_consecutive = k_consecutive  # Number of consecutive reads required to support a decision
        self.consecutive_reads = defaultdict(list)  # Store last reads info for each line/chromosome combination
        self.current_haplotype_state = defaultdict(lambda: {'haplotype': None, 'count': 0})
        self.previous_phase_set = None


    def calculate_score(self, read):
        hp1_count, hp2_count = 0, 0
        haplotype_blocks = []
        variant_string = []
        haplotype_string = []
        current_haplotype = None
        block_start = None
        block_end = None
        block_size = 0
        ref_pos = read.reference_start
        query_pos = 0
        ps_list = []

        try:
            line = read.get_tag("LN")
        except KeyError:
            line = "none"

        for cigar_op, length in read.cigartuples:
            if cigar_op in [0, 7, 8]:  # Match, =, X (regular aligned bases)
                for i in range(length):
                    pos = ref_pos + i
                    if read.reference_name in self.variants and pos in self.variants[read.reference_name]:
                        if read.query_qualities[query_pos + i-1] < self.min_variant_quality:  # SNP-specific quality check
                            variant_string.append(
                                f"{pos}:{read.query_sequence[query_pos + i - 1]}:{read.query_qualities[query_pos + i - 1]}:quality")
                            continue

                        var = self.variants[read.reference_name][pos]
                        if len(var['ref']) == 1 and len(var['alt']) == 1:  # SNP case
                            # Process SNPs
                            hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size, variant_str, haplotype_str, phase_set = self._process_snp(
                                read, query_pos + i - 1, var, hp1_count, hp2_count, current_haplotype, block_start,
                                block_end, pos, haplotype_blocks, block_size
                            )
                            ps_list.append(phase_set)
                            variant_string.append(variant_str)
                            haplotype_string.append(haplotype_str)

                query_pos += length
                ref_pos += length

            elif cigar_op == 1 or cigar_op == 2:  # Insertions or Deletions
                pos = ref_pos if cigar_op == 2 else ref_pos   # For deletions, pos should be at ref_pos; for insertions, it's at the prior position.
                if read.reference_name in self.variants and pos in self.variants[read.reference_name]:
                    var = self.variants[read.reference_name][pos]

                    # Process Indels (either insertion or deletion based on cigar_op)
                    hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size, variant_str, haplotype_str, phase_set = self._process_indel(
                        read, query_pos , var, hp1_count, hp2_count, current_haplotype, block_start, block_end, pos,
                        haplotype_blocks, block_size, cigar_op, length
                    )
                    ps_list.append(phase_set)
                    variant_string.append(variant_str)
                    haplotype_string.append(haplotype_str)

                # Advance positions after indels (insertion affects query position, deletion affects reference position)
                if cigar_op == 1:  # Insertion
                    query_pos += length
                elif cigar_op == 2:  # Deletion
                    ref_pos += length

            elif cigar_op == 4:  # Soft clipping
                query_pos += length

        if current_haplotype is not None and block_size >= self.min_supporting_variants:
            haplotype_blocks.append((block_start, block_end, current_haplotype, block_size))

        # Calculate the dominant haplotype for the current read
        dominant_haplotype, dominant_fraction = self.calculate_dominant_haplotype("".join(haplotype_string))
        switch_score = self.calculate_phase_switches("".join(haplotype_string))
        con_score = self.calculate_consecutiveness_score("".join(haplotype_string))

        dominant_phase_set = max(ps_list, key=ps_list.count) if ps_list else None

        # Get chromosome and line for storing in the dictionary
        chrom_line_key = (read.reference_name, line)

        # Check if the read is recombinant
        (is_recombinant, breakpoints, previous_read_id,
         scored_read_id,scored_read_start,scored_read_haplotype_blocks,
         scored_read_variant_string,scored_read_haplotype_string,
         scored_read_switch_score,scored_read_con_score
         )= (False, [], None, None,None,[],[],[],None,None ) if not haplotype_blocks else self.is_recombinant(
            chrom_line_key, haplotype_blocks, switch_score, con_score, dominant_haplotype, dominant_fraction,
            dominant_phase_set, read.query_name, line, read.query_name, read.reference_start, variant_string,
            haplotype_string
        )

        # Store the current read's info if it has a valid dominant haplotype
        if dominant_haplotype is not None and dominant_fraction >= self.dominant_fraction_threshold and line != "none":
            self.previous_read_info[chrom_line_key] = {
                'dominant_haplotype': dominant_haplotype,
                'haplotype_blocks': haplotype_blocks,
                'read_id': read.query_name,
                'dominant_phase_set': dominant_phase_set
            }

        return {
            'read_id': scored_read_id,
            'chrom': read.reference_name,
            'start_pos': scored_read_start,
            'is_recombinant': is_recombinant,
            'previous_read_id': previous_read_id,  # Include the previous read ID in the output
            'haplotype_blocks': scored_read_haplotype_blocks,
            'variant_string': ";".join(scored_read_variant_string),
            'haplotype_string': "".join(scored_read_haplotype_string),
            'breakpoints': breakpoints,
            'line': line,
            'switch_score': scored_read_switch_score,
            'con_score': scored_read_con_score
        }

    def is_recombinant(self, chrom_line_key, haplotype_blocks, switch_score, con_score,
                       dominant_haplotype, dominant_fraction, dominant_phase_set, read_name, line,
                       read_id, read_start, variant_string, haplotype_string):
        """Determines if the current read is recombinant, based on single-read and multi-read logic."""
        # Initialize recombination status
        is_recombinant = False
        previous_read_id = None
        breakpoints = []

        # Check for single-read recombination
        (single_read_recombinant, single_read_breakpoints,
            previous_read_id, read_id, read_start, haplotype_blocks,
         variant_string, haplotype_string, switch_score, con_score
         )= self.is_single_read_recombinant(haplotype_blocks, switch_score, con_score, previous_read_id, read_id, read_start, variant_string, haplotype_string)
        if single_read_recombinant:
            is_recombinant = True
            breakpoints = single_read_breakpoints

        # Check for multi-read recombination if no single-read recombination occurred
        if not is_recombinant and dominant_haplotype is not None and line != "none" and dominant_fraction >= self.dominant_fraction_threshold:
            (multi_read_recombinant, multi_read_breakpoints,
            previous_read_id, read_id, read_start, haplotype_blocks,
            variant_string, haplotype_string, switch_score, con_score
            )= self.is_multi_read_recombinant(
                chrom_line_key, dominant_haplotype, haplotype_blocks, read_name, dominant_phase_set, read_start, variant_string,
                haplotype_string, switch_score, con_score
            )
            if multi_read_recombinant:
                is_recombinant = True
                breakpoints = multi_read_breakpoints

        return(is_recombinant, breakpoints,
        previous_read_id,
        read_id,
        read_start,
        haplotype_blocks,
        variant_string,
        haplotype_string,
        switch_score,
        con_score)

    def is_single_read_recombinant(self, haplotype_blocks, switch_score, con_score,
                                   previous_read_id, read_id, read_start, variant_string, haplotype_string):
        """Checks if the recombination event is within the single read."""
        is_recombinant = False
        breakpoints = []

        for i in range(len(haplotype_blocks) - 1):
            if (haplotype_blocks[i][3] >= self.min_supporting_variants and
                    haplotype_blocks[i + 1][3] >= self.min_supporting_variants and
                    haplotype_blocks[i][2] != haplotype_blocks[i + 1][2]) and (
                    switch_score < self.switch_threshold and con_score > self.continuity_threshold):
                is_recombinant = True
                end_pos1 = haplotype_blocks[i][1]
                start_pos2 = haplotype_blocks[i + 1][0]
                breakpoints.append(f"{end_pos1}-{start_pos2}")

        return (is_recombinant, breakpoints,
                previous_read_id,
                read_id,
                read_start,
                haplotype_blocks,
                variant_string,
                haplotype_string,
                switch_score,
                con_score)

    def is_multi_read_recombinant(self, chrom_line_key, current_haplotype, current_haplotype_blocks,
                                  read_name, dominant_phase_set, read_start, variant_string,
                                  haplotype_string, switch_score, con_score):
        """Checks if the recombination event involves consecutive reads before and after a switch and calculates the correct breakpoints."""
        is_recombinant = False
        previous_read_id = None
        breakpoints = []
        end = None
        start = None

        # Clean up `consecutive_reads` if the new read has a different `dominant_phase_set`
        if (chrom_line_key in self.consecutive_reads and self.consecutive_reads[chrom_line_key] and
                self.consecutive_reads[chrom_line_key][-1]['dominant_phase_set'] != dominant_phase_set):
            self.consecutive_reads[chrom_line_key] = []

        # Extract the last variant position from the appropriate haplotype block
        for block in reversed(current_haplotype_blocks):
            if block[2][2] == current_haplotype:  # Ensure we get the last block corresponding upstream haplotype
                end = block[1]  # Take the 'end' position of this block
                break

        # Extract the first variant position from the appropriate haplotype block
        for block in current_haplotype_blocks:
            if block[2][2] == current_haplotype:  # Ensure we get the first block corresponding downstream haplotype
                start = block[0]  # Take the 'start' position of this block
                break


        # Store the current read in the list of consecutive reads
        self.consecutive_reads[chrom_line_key].append({
            'dominant_haplotype': current_haplotype,
            'dominant_phase_set': dominant_phase_set,
            'haplotype_blocks': current_haplotype_blocks,
            'read_id': read_name,
            'start': start,
            'end': end,
            'read_start': read_start,
            'read_variant_string': variant_string,
            'read_haplotype_string': haplotype_string,
            'read_switch_score': switch_score,
            'read_con_score': con_score
        })

        # Ensure we are only looking at the last `2 * k_consecutive` reads
        if len(self.consecutive_reads[chrom_line_key]) > 2 * self.k_consecutive:
            self.consecutive_reads[chrom_line_key].pop(0)

        # Check if we have enough consecutive reads (before and after) to analyze
        if len(self.consecutive_reads[chrom_line_key]) < 2 * self.k_consecutive:
            return (is_recombinant, breakpoints, previous_read_id,
                    read_name, read_start, current_haplotype_blocks, variant_string,
                    haplotype_string, switch_score, con_score
                    )


        scored_read = self.consecutive_reads[chrom_line_key][self.k_consecutive]
        scored_read_id = scored_read['read_id']
        scored_read_start = scored_read['read_start']
        scored_read_haplotype_blocks = scored_read['haplotype_blocks']
        scored_read_variant_string = scored_read['read_variant_string']
        scored_read_haplotype_string = scored_read['read_haplotype_string']
        scored_read_switch_score = scored_read['read_switch_score']
        scored_read_con_score = scored_read['read_con_score']

        # Extract the reads before, during, and after the potential crossover event
        reads_before = self.consecutive_reads[chrom_line_key][:self.k_consecutive-1]
        reads_after = self.consecutive_reads[chrom_line_key][self.k_consecutive:]

        # Check if there is a stable haplotype before and after the middle read
        stable_before = all(
            read['dominant_haplotype'] == reads_before[0]['dominant_haplotype'] for read in reads_before)
        stable_after = all(read['dominant_haplotype'] == reads_after[0]['dominant_haplotype'] for read in reads_after)

        if stable_before and stable_after and reads_before[0]['dominant_haplotype'] != reads_after[0][
            'dominant_haplotype'] and reads_before[0]['dominant_phase_set'] == reads_after[0]['dominant_phase_set']:
            # If the haplotype switches and remains stable both before and after the switch
            is_recombinant = True


            # Calculate the breakpoints
            breakpoints = self._calculate_breakpoints(reads_before, reads_after)
            previous_read_id = breakpoints[2]
            breakpoints = [f"{breakpoints[0]}-{breakpoints[1]}"]



        return (is_recombinant, breakpoints,
        previous_read_id,
        scored_read_id,
        scored_read_start,
        scored_read_haplotype_blocks,
        scored_read_variant_string,
        scored_read_haplotype_string,
        scored_read_switch_score,
        scored_read_con_score)

    def _calculate_breakpoints(self, reads_before, reads_after):
        """
        Calculates the breakpoints between haplotype 1 and haplotype 2.

        Breakpoint should be between the last variant in the last read of haplotype 1
        and the first variant in the first read of haplotype 2.
        """
        # Handle ties by read ID or other features if necessary
        last_read_before = max(reads_before, key=lambda x: (x['end'], x['read_id']))
        last_haplotype1_position, last_read_id = last_read_before['end'], last_read_before[
            'read_id']  # End of the last read covering haplotype 1

        first_read_after = min(reads_after, key=lambda x: (x['start'], x['read_id']))
        first_haplotype2_position = first_read_after['start']  # Start of the first read covering haplotype 2

        # The breakpoint is between the last haplotype 1 position and the first haplotype 2 position
        return [last_haplotype1_position, first_haplotype2_position, last_read_id]



    @staticmethod
    def calculate_dominant_haplotype(haplotype_string):
        """Calculates the dominant haplotype and its fraction."""
        count_1 = haplotype_string.count('1')
        count_2 = haplotype_string.count('2')
        total = count_1 + count_2

        if count_1 < 3 and count_2 < 3:
            return None, 0

        if total == 0:
            return None, 0

        dominant_haplotype = '1' if count_1 > count_2 else '2'
        dominant_fraction = max(count_1, count_2) / total

        return dominant_haplotype, dominant_fraction

    @staticmethod
    def calculate_phase_switches(haplotype_string):
        """Calculates the phase switches in the haplotype string."""
        switches = sum(1 for i in range(1, len(haplotype_string)) if
                       haplotype_string[i] != haplotype_string[i - 1] and haplotype_string[i] in '12')
        return switches / len(haplotype_string) if haplotype_string else 0

    @staticmethod

    ## It's a shitty metric, think about squaring the score to give more weight to longer consecutive blocks
    ## And normalizing by the length of the haplotype string or number of blocks
    def calculate_consecutiveness_score(haplotype_string):
        def score_for_haplotype(hap):
            score, consecutive, total_count = 0, 0, haplotype_string.count(hap)
            for char in haplotype_string:
                if char == hap:
                    consecutive += 1
                    score += consecutive
                elif char == 'x':
                    continue
                else:
                    consecutive = 0
            return score / total_count if total_count else 0

        score_1 = score_for_haplotype('1')
        score_2 = score_for_haplotype('2')
        average_score = (score_1 + score_2) / 2
        return average_score / len(haplotype_string) if haplotype_string else 0

    def _process_snp(self, read, pos, var, hp1_count, hp2_count, current_haplotype, block_start, block_end, ref_pos,
                     haplotype_blocks, block_size):
        """
        Processes a single SNP (Single Nucleotide Polymorphism) and updates haplotype counts and blocks accordingly.

        Arguments:
        - read: The read being processed
        - pos: Position in the read's query sequence
        - var: The variant information (from VCF) at this position
        - hp1_count: Count of hp1 variants in the read
        - hp2_count: Count of hp2 variants in the read
        - current_haplotype: The current haplotype being followed in this block (hp1 or hp2)
        - block_start: Start of the current haplotype block
        - block_end: End of the current haplotype block
        - ref_pos: The reference position being processed
        - haplotype_blocks: List of haplotype blocks identified so far
        - block_size: Number of variants in the current haplotype block

        Returns:
        Updated counts and block information after processing this SNP.
        """
        read_base = read.query_sequence[pos]
        variant_str = f"{ref_pos}:{read.query_sequence[pos]}:{'hp1' if read_base == var['hp1'] else 'hp2'}:{read.query_qualities[pos]}:"

        # Phase set must match to continue processing the current haplotype
        current_phase_set = var['ps']

        if read_base == var['hp1']:
            haplotype_str = '1'
        elif read_base == var['hp2']:
            haplotype_str = '2'
        else:
            haplotype_str = 'x'

        # Check for haplotype switch or phase set switch
        if read_base == var['hp1']:
            hp1_count += 1
            if current_haplotype != 'hp1' or (block_size > 0 and current_phase_set != self.previous_phase_set): ## fix previous_phase_set.
                if current_haplotype is not None and block_size >= self.min_supporting_variants:
                    haplotype_blocks.append((block_start, block_end, current_haplotype, block_size))
                current_haplotype = 'hp1'
                block_start = ref_pos
                block_size = 0
            block_end = ref_pos
            block_size += 1
        elif read_base == var['hp2']:
            hp2_count += 1
            if current_haplotype != 'hp2' or (block_size > 0 and current_phase_set != self.previous_phase_set):
                if current_haplotype is not None and block_size >= self.min_supporting_variants:
                    haplotype_blocks.append((block_start, block_end, current_haplotype, block_size))
                current_haplotype = 'hp2'
                block_start = ref_pos
                block_size = 0
            block_end = ref_pos
            block_size += 1

        # Update the previous phase set to the current one
        self.previous_phase_set = current_phase_set

        return hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size, variant_str, haplotype_str, \
        var['ps']

    def _process_indel(self, read, query_pos, variant, hp1_count, hp2_count, current_haplotype, block_start, block_end,
                       pos, haplotype_blocks, block_size, cigar_op,length):
        """
        Processes an indel (insertion or deletion) and handles sequence comparison, quality checks, haplotype phase, and block management.

        Args:
            - cigar_op: Indicates whether it's an insertion (1) or deletion (2)
        """
        ref_seq = variant['ref']
        alt_seq = variant['alt']
        tolerance_percentage = 0.1  # Configurable tolerance

        # Handle insertions
        if cigar_op == 1:  # Insertion
            read_seq = read.query_sequence[query_pos-1:query_pos-1 + len(alt_seq)]
            inserted_qualities = read.query_qualities[query_pos-1:query_pos-1 + len(alt_seq)]
            mean_quality = sum(inserted_qualities) / len(inserted_qualities)

            if mean_quality < self.min_variant_quality:
                return hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size, f"{pos}:+{alt_seq}:low_quality", "", None

            # Sequence comparison with tolerance using Levenshtein (for insertions)
            levenshtein_dist = levenshtein_distance(read_seq, alt_seq)
            max_mismatches = int(len(alt_seq) * tolerance_percentage)

            if levenshtein_dist > max_mismatches:
                #print(variant, file=sys.stderr)
                #print(f"{read.query_name} {pos} {read_seq} {alt_seq} {levenshtein_dist} {max_mismatches}", file=sys.stderr)
                return hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size, f"{pos}:+{alt_seq}:mismatch", "", None

        # Handle deletions
        elif cigar_op == 2:  # Deletion
            deletion_length = len(ref_seq)-1

            # Calculate allowed wiggle room based on tolerance_percentage
            tolerance = int(deletion_length * tolerance_percentage)
            min_length = deletion_length - tolerance
            max_length = deletion_length + tolerance

            # If the CIGAR length is outside the allowed range, report a length mismatch
            if not (min_length <= length <= max_length):  # `length` comes from the CIGAR tuple
                #print(variant, file=sys.stderr)
                #print (f"{read.query_name} {pos} {deletion_length} {length} {min_length} {max_length}", file=sys.stderr)
                return hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size, f"{pos}:-{ref_seq}:length_mismatch", "", None

        # Determine the haplotype from the variant (not the read!)
        if alt_seq == variant['hp1']:
            haplotype_str = '1'
            hp1_count += 1
            if current_haplotype != 'hp1' or (block_size > 0 and variant['ps'] != self.previous_phase_set):
                if current_haplotype is not None and block_size >= self.min_supporting_variants:
                    haplotype_blocks.append((block_start, block_end, current_haplotype, block_size))
                current_haplotype = 'hp1'
                block_start = pos
                block_size = 0
            block_end = pos
            block_size += 1
        elif alt_seq == variant['hp2']:
            haplotype_str = '2'
            hp2_count += 1
            if current_haplotype != 'hp2' or (block_size > 0 and variant['ps'] != self.previous_phase_set):
                if current_haplotype is not None and block_size >= self.min_supporting_variants:
                    haplotype_blocks.append((block_start, block_end, current_haplotype, block_size))
                current_haplotype = 'hp2'
                block_start = pos
                block_size = 0
            block_end = pos
            block_size += 1
        else:
            haplotype_str = 'x'  # Unphased or unrecognized haplotype

        # Update the previous phase set
        self.previous_phase_set = variant['ps']

        # Construct variant string for reporting
        variant_str = f"{pos}:{'+'+alt_seq if cigar_op == 1 else '-'+ref_seq}:{'hp'+haplotype_str}:{variant['ps']}"

        return hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size, variant_str, haplotype_str, \
            variant['ps']


def defaultdict_of_dict():
    return defaultdict(dict)

def levenshtein_distance(seq1, seq2):
    """
    Calculates the Levenshtein distance between two sequences (seq1 and seq2).
    This metric allows insertion, deletion, or substitution with a cost of 1 for each.
    """
    len1, len2 = len(seq1), len(seq2)
    matrix = [[0] * (len2 + 1) for _ in range(len1 + 1)]

    for i in range(len1 + 1):
        matrix[i][0] = i
    for j in range(len2 + 1):
        matrix[0][j] = j

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if seq1[i - 1] == seq2[j - 1]:
                cost = 0
            else:
                cost = 1
            matrix[i][j] = min(matrix[i - 1][j] + 1,      # Deletion
                               matrix[i][j - 1] + 1,      # Insertion
                               matrix[i - 1][j - 1] + cost)  # Substitution

    return matrix[len1][len2]


def parse_vcf(vcf_file):
    variants = defaultdict(defaultdict_of_dict)
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            chrom, pos, _, ref, alt, _, _, _, fmt, sample = line.strip().split('\t')
            pos = int(pos)
            gt = sample.split(':')[0]
            ps = sample.split(':')[-1]
            hp1, hp2 = gt.split('|')
            if chrom not in variants:
                variants[chrom] = {}
            var = {
                'ref': ref,
                'alt': alt,
                'hp1': ref if hp1 == '0' else alt,
                'hp2': ref if hp2 == '0' else alt,
                'ps': int(ps)
            }
            variants[chrom][pos] = var
    return variants


def main():
    parser = argparse.ArgumentParser(description="Detect recombination events in reads.")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file with phased variants.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file with aligned reads.")
    parser.add_argument("-m", "--min_supporting_variants", type=int, default=5,
                        help="Minimum number of supporting variants to consider a phase switch.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of worker processes to use.")
    parser.add_argument("-o", "--output", help="Output file to write results. Writes to stdout if not specified.")
    parser.add_argument("-c", "--continuity_threshold", type=float, default=0.15,
                        help="Threshold for continuity score to determine recombination.")
    parser.add_argument("-s", "--switch_threshold", type=float, default=0.1,
                        help="Threshold for switch score to determine recombination.")
    parser.add_argument("--dominant_fraction_threshold", type=float, default=0.9,
                        help="Fraction threshold for determining dominant haplotype.")
    parser.add_argument("--min_variant_quality", type=int, default=20, help="Minimum quality score for variants.")
    parser.add_argument("--min_mapping_quality", type=int, default=50, help="Minimum mapping quality for reads.")
    parser.add_argument("--min_read_length", type=int, default=1000, help="Minimum read length.")
    args = parser.parse_args()

    variants = parse_vcf(args.vcf)
    results = process_bam(args.bam, variants, args.min_supporting_variants, args.continuity_threshold,
                          args.switch_threshold, args.dominant_fraction_threshold, args.min_mapping_quality,
                          args.min_read_length, args.min_variant_quality, num_workers=args.threads)

    output_file = sys.stdout if args.output is None else open(args.output, 'w')

    # Updated header to include previous_read_id
    header = "ReadID\tLine\tChromosome\tStartPos\tIsRecombinant\tPreviousReadID\tHaplotypeBlocks\tVariantString\tHaplotypeString\tBreakpoints\tSwitchScore\tConScore\n"
    output_file.write(header)

    for result in results:
        haplotype_blocks_str = ";".join([f"{block[0]}-{block[1]}:{block[2]}:{block[3]}" for block in result['haplotype_blocks']])
        breakpoints_str = ";".join(map(str, result['breakpoints']))


        line = (
            f"{result['read_id']}\t"
            f"{result['line']}\t"
            f"{result['chrom']}\t"
            f"{result['start_pos']}\t"
            f"{result['is_recombinant']}\t"
            f"{result['previous_read_id'] if result['is_recombinant'] and result['previous_read_id'] else 'NA'}\t"  # Handle previous_read_id
            f"{haplotype_blocks_str}\t"
            f"{result['variant_string']}\t"
            f"{result['haplotype_string']}\t"
            f"{breakpoints_str}\t"
            f"{result['switch_score']}\t"
            f"{result['con_score']}\n"
        )
        output_file.write(line)

    if args.output is not None:
        output_file.close()


if __name__ == "__main__":
    main()

