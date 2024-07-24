#! /usr/bin/env python

import pysam    # Import the pysam module
import argparse     # Import the argparse module
from collections import defaultdict     # Import the defaultdict class from the collections module
import sys      # Import the sys module

class HaplotypeAnalyzer:
    """
    Analyzes haplotypes in reads and calculates haplotype scores.
    """
    def __init__(self, variants, min_supporting_variants):
        """
        Initializes the HaplotypeAnalyzer with given variants and minimum supporting variants.
        """
        self.variants = variants
        self.min_supporting_variants = min_supporting_variants

    def calculate_score(self, read):
        """
        Calculates the haplotype score for a given read.
        """
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

        # Iterate over the CIGAR operations
        for cigar_op, length in read.cigartuples:
            if cigar_op in [0, 7, 8]:  # Match, =, X
                for i in range(length):
                    pos = ref_pos + i
                    if pos in self.variants[read.reference_name]:
                        if read.query_qualities[query_pos + i - 1] < 20:
                            variant_string.append(f"{pos}:{read.query_sequence[query_pos + i - 1]}:{read.query_qualities[query_pos + i - 1]}:quality")
                            continue
                        var = self.variants[read.reference_name][pos]
                        if len(var['ref']) == 1 and len(var['alt']) == 1:  # SNP
                            hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size, variant_str, haplotype_str = self._process_snp(
                                read, query_pos + i - 1, var, hp1_count, hp2_count, current_haplotype, block_start, block_end, pos, haplotype_blocks, block_size
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

        # Finalize the last haplotype block if exists
        if current_haplotype is not None and block_size >= self.min_supporting_variants:
            haplotype_blocks.append((block_start, block_end, current_haplotype, block_size))

        is_recombinant, breakpoints = self._finalize_score(haplotype_blocks)

        switch_score = self.calculate_phase_switches("".join(haplotype_string))
        con_score = self.calculate_consecutiveness_score("".join(haplotype_string))

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

    def _process_snp(self, read, pos, var, hp1_count, hp2_count, current_haplotype, block_start, block_end, ref_pos, haplotype_blocks, block_size):
        """
        Processes a single SNP variant.
        """
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



        return (hp1_count, hp2_count, current_haplotype, block_start, block_end, block_size,
                variant_str, haplotype_str)

    def _finalize_score(self, haplotype_blocks):
        """
        Finalizes the haplotype score by determining if the read is recombinant.
        Collects breakpoints between consecutive off-phase blocks for recombinant reads.
        """
        is_recombinant = False
        breakpoints = []
        
        for i in range(len(haplotype_blocks) - 1):
            if (haplotype_blocks[i][3] >= self.min_supporting_variants and 
                haplotype_blocks[i + 1][3] >= self.min_supporting_variants and
                haplotype_blocks[i][2] != haplotype_blocks[i + 1][2]):
                
                is_recombinant = True
                
                # Determine innermost variants between consecutive off-phase blocks
                end_pos1 = haplotype_blocks[i][1]
                start_pos2 = haplotype_blocks[i + 1][0]
                
                breakpoints.append(f"{end_pos1}-{start_pos2}")
        
        return is_recombinant, breakpoints

    def calculate_phase_switches(self, haplotype_string):
        switches = sum(1 for i in range(1, len(haplotype_string)) if
                       haplotype_string[i] != haplotype_string[i-1] and haplotype_string[i] in '12')
        return switches / len(haplotype_string) if haplotype_string else 0

    def calculate_consecutiveness_score(self, haplotype_string):
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

def parse_vcf(vcf_file):
    """
    Parses the VCF file and extracts variants.
    """
    variants = defaultdict(lambda: defaultdict(dict))
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            chrom, pos, _, ref, alt, _, _, _, fmt, sample = line.strip().split('\t')
            pos = int(pos)
            gt = sample.split(':')[0]
            ps = sample.split(':')[-1]
            hp1, hp2 = gt.split('|')
            var = {
                'ref': ref,
                'alt': alt,
                'hp1': ref if hp1 == '0' else alt,
                'hp2': ref if hp2 == '0' else alt,
                'ps': int(ps)
            }
            variants[chrom][pos] = var
    return variants

def process_bam(bam_file, variants, min_supporting_variants):
    """
    Processes the BAM file and calculates haplotype scores for each read.
    """
    analyzer = HaplotypeAnalyzer(variants, min_supporting_variants)
    results = []

    bam = pysam.AlignmentFile(bam_file, "rb")
    for read in bam.fetch():
        if read.is_secondary or read.is_supplementary:
            continue

        result = analyzer.calculate_score(read)
        results.append(result)

    return results

def main():
    """
    Main function to parse arguments, process VCF and BAM files, and output results.
    """
    parser = argparse.ArgumentParser(description="Detect recombination events in reads.")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file with phased variants.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file with aligned reads.")
    parser.add_argument("-m", "--min_supporting_variants", type=int, default=2, help="Minimum number of supporting variants to consider a phase switch.")
    parser.add_argument("-o", "--output", help="Output file to write results. Writes to stdout if not specified.")

    args = parser.parse_args()

    # Capture the command used to run the tool
    command_used = ' '.join(sys.argv)

    variants = parse_vcf(args.vcf)
    results = process_bam(args.bam, variants, args.min_supporting_variants)

    output_file = sys.stdout if args.output is None else open(args.output, 'w')

    # Write the command used as a header in the output file
    output_file.write(f"# Command: {command_used}\n")

    header = "ReadID\tLine\tChromosome\tStartPos\tIsRecombinant\tHaplotypeBlocks\tVariantString\tHaplotypeString\tBreakpoints\tSwitchScore\tConScore\n"
    output_file.write(header)
    for result in results:
        is_recombinant, breakpoints = result['is_recombinant'], result.get('breakpoints', [])

        output_file.write(
            f"{result['read_id']}\t"
            f"{result['line']}\t"
            f"{result['chrom']}\t"
            f"{result['start_pos']}\t"
            f"{is_recombinant}\t"
            f"{result['haplotype_blocks']}\t"
            f"{result['variant_string']}\t"
            f"{result['haplotype_string']}\t"
            f"{';'.join(breakpoints)}\t"
            f"{result['switch_score']}\t"
            f"{result['con_score']}\n"
        )

    if args.output is not None:
        output_file.close()

if __name__ == "__main__":
    main()
