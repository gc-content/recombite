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
            end = start + chunk_size
            if end > length:
                end = length
            # Adjust the start of the next range to avoid overlap
            if end == length:
                chromosome_ranges.append((chromosome, start, end))
            else:
                chromosome_ranges.append((chromosome, start, end - 1))

    bam.close()
    return chromosome_ranges


def process_bam_chunk(bam_file, variants, min_supporting_variants, continuity_threshold, switch_threshold, chunk):
    chromosome, start, end = chunk
    analyzer = HaplotypeAnalyzer(variants, min_supporting_variants, continuity_threshold, switch_threshold)
    results = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chromosome, start, end):
            if read.is_secondary or read.is_supplementary:
                continue
            result = analyzer.calculate_score(read)
            results.append(result)

    return results


def process_bam(bam_file, variants, min_supporting_variants, continuity_threshold, switch_threshold, num_workers=4):
    chromosome_ranges = get_chromosome_ranges(bam_file, num_workers)

    with multiprocessing.Pool(processes=num_workers) as pool:
        results = pool.starmap(
            process_bam_chunk,
            [(bam_file, variants, min_supporting_variants, continuity_threshold, switch_threshold, chunk) for chunk in chromosome_ranges]
        )

    # Flatten the list of results
    return [item for sublist in results for item in sublist]


class HaplotypeAnalyzer:
    def __init__(self, variants, min_supporting_variants, continuity_threshold, switch_threshold):
        self.variants = variants
        self.min_supporting_variants = min_supporting_variants
        self.continuity_threshold = continuity_threshold
        self.switch_threshold = switch_threshold

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

        try:
            line = read.get_tag("LN")
        except KeyError:
            line = "none"

        for cigar_op, length in read.cigartuples:
            if cigar_op in [0, 7, 8]:  # Match, =, X
                for i in range(length):
                    pos = ref_pos + i
                    if read.reference_name in self.variants and pos in self.variants[read.reference_name]:
                        if read.query_qualities[query_pos + i - 1] < 20:
                            variant_string.append(
                                f"{pos}:{read.query_sequence[query_pos + i - 1]}:{read.query_qualities[query_pos + i - 1]}:quality")
                            continue
                        var = self.variants[read.reference_name][pos]
                        if len(var['ref']) == 1 and len(var['alt']) == 1:  # SNP
                            hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size, variant_str, haplotype_str = self._process_snp(
                                read, query_pos + i - 1, var, hp1_count, hp2_count, current_haplotype, block_start,
                                block_end, pos, haplotype_blocks, block_size
                            )
                            variant_string.append(variant_str)
                            haplotype_string.append(haplotype_str)
                query_pos += length
                ref_pos += length

            elif cigar_op == 1:  # Insertion
                query_pos += length

            elif cigar_op == 2:  # Deletion
                ref_pos += length

            elif cigar_op == 4:  # Soft clipping
                query_pos += length

        if current_haplotype is not None and block_size >= self.min_supporting_variants:
            haplotype_blocks.append((block_start, block_end, current_haplotype, block_size))


        switch_score = self.calculate_phase_switches("".join(haplotype_string))
        con_score = self.calculate_consecutiveness_score("".join(haplotype_string))

        is_recombinant, breakpoints = self._finalize_score(haplotype_blocks, switch_score, con_score)

        return {
            'read_id': read.query_name,
            'chrom': read.reference_name,
            'start_pos': read.reference_start,
            'is_recombinant': is_recombinant,
            'haplotype_blocks': haplotype_blocks,
            'variant_string': ";".join(variant_string),
            'haplotype_string': "".join(haplotype_string),
            'breakpoints': breakpoints,
            'line': line,
            'switch_score': switch_score,
            'con_score': con_score
        }

    def _process_snp(self, read, pos, var, hp1_count, hp2_count, current_haplotype, block_start, block_end, ref_pos,
                     haplotype_blocks, block_size):
        read_base = read.query_sequence[pos]
        variant_str = f"{ref_pos}:{read.query_sequence[pos]}:{'hp1' if read_base == var['hp1'] else 'hp2'}:{read.query_qualities[pos]}:"

        if read_base == var['hp1']:
            haplotype_str = '1'
        elif read_base == var['hp2']:
            haplotype_str = '2'
        else:
            haplotype_str = 'x'

        if read_base == var['hp1']:
            hp1_count += 1
            if current_haplotype != 'hp1':
                if current_haplotype is not None and block_size >= self.min_supporting_variants:
                    haplotype_blocks.append((block_start, block_end, current_haplotype, block_size))
                current_haplotype = 'hp1'
                block_start = ref_pos
                block_size = 0
            block_end = ref_pos
            block_size += 1
        elif read_base == var['hp2']:
            hp2_count += 1
            if current_haplotype != 'hp2':
                if current_haplotype is not None and block_size >= self.min_supporting_variants:
                    haplotype_blocks.append((block_start, block_end, current_haplotype, block_size))
                current_haplotype = 'hp2'
                block_start = ref_pos
                block_size = 0
            block_end = ref_pos
            block_size += 1

        return hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size, variant_str, haplotype_str

    def _finalize_score(self, haplotype_blocks, switch_score, con_score):
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

        return is_recombinant, breakpoints

    @staticmethod
    def calculate_phase_switches(haplotype_string):
        switches = sum(1 for i in range(1, len(haplotype_string)) if
                       haplotype_string[i] != haplotype_string[i - 1] and haplotype_string[i] in '12')
        return switches / len(haplotype_string) if haplotype_string else 0

    @staticmethod
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


def defaultdict_of_dict():
    return defaultdict(dict)


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
    parser.add_argument("-m", "--min_supporting_variants", type=int, default=2,
                        help="Minimum number of supporting variants to consider a phase switch.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of worker processes to use.")
    parser.add_argument("-o", "--output", help="Output file to write results. Writes to stdout if not specified.")
    parser.add_argument("-c", "--continuity_threshold", type=float, default=0.2,
                        help="Threshold for continuity score to determine recombination.")
    parser.add_argument("-s", "--switch_threshold", type=float, default=0.01,
                        help="Threshold for switch score to determine recombination.")
    args = parser.parse_args()

    variants = parse_vcf(args.vcf)
    results = process_bam(args.bam, variants, args.min_supporting_variants, args.continuity_threshold,
                          args.switch_threshold, num_workers=args.threads)

    output_file = sys.stdout if args.output is None else open(args.output, 'w')

    header = "ReadID\tLine\tChromosome\tStartPos\tIsRecombinant\tHaplotypeBlocks\tVariantString\tHaplotypeString\tBreakpoints\tSwitchScore\tConScore\n"
    output_file.write(header)

    for result in results:
        haplotype_blocks_str = ";".join([f"{block[0]}-{block[1]}:{block[2]}:{block[3]}" for block in result['haplotype_blocks']])
        breakpoints_str = ";".join(result['breakpoints'])

        line = (
            f"{result['read_id']}\t"
            f"{result['line']}\t"
            f"{result['chrom']}\t"
            f"{result['start_pos']}\t"
            f"{result['is_recombinant']}\t"
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
