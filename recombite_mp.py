import pysam
import argparse
import sys
import multiprocessing

def get_chromosome_ranges(bam_file, num_chunks):
    """
    Divide the BAM file's reference chromosomes into ranges for parallel processing.
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

def process_bam_chunk(bam_file, variants, min_supporting_variants, chunk):
    """
    Process a specific chunk of the BAM file.
    """
    chromosome, start, end = chunk
    analyzer = HaplotypeAnalyzer(variants, min_supporting_variants)
    results = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chromosome, start, end):
            if read.is_secondary or read.is_supplementary:
                continue
            result = analyzer.calculate_score(read)
            results.append(result)

    return results

def process_bam(bam_file, variants, min_supporting_variants, num_workers=4):
    """
    Main function to process the BAM file with multiprocessing.
    """
    chromosome_ranges = get_chromosome_ranges(bam_file, num_workers)
    
    with multiprocessing.Pool(processes=num_workers) as pool:
        results = pool.starmap(
            process_bam_chunk,
            [(bam_file, variants, min_supporting_variants, chunk) for chunk in chromosome_ranges]
        )
    
    # Flatten the list of results
    return [item for sublist in results for item in sublist]

class HaplotypeAnalyzer:
    def __init__(self, variants, min_supporting_variants):
        self.variants = variants
        self.min_supporting_variants = min_supporting_variants

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

        if current_haplotype is not None and block_size >= self.min_supporting_variants:
            haplotype_blocks.append((block_start, block_end, current_haplotype, block_size))

        is_recombinant, breakpoints = self._finalize_score(haplotype_blocks)

        return {
            'read_id': read.query_name,
            'chrom': read.reference_name,
            'start_pos': read.reference_start,
            'is_recombinant': is_recombinant,
            'haplotype_blocks': haplotype_blocks,
            'variant_string': ";".join(variant_string),
            'haplotype_string': "".join(haplotype_string),
            'breakpoints': breakpoints,
            'line': line
        }

    def _process_snp(self, read, pos, var, hp1_count, hp2_count, current_haplotype, block_start, block_end, ref_pos, haplotype_blocks, block_size):
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

    def _finalize_score(self, haplotype_blocks):
        is_recombinant = False
        breakpoints = []
        
        for i in range(len(haplotype_blocks) - 1):
            if (haplotype_blocks[i][3] >= self.min_supporting_variants and 
                haplotype_blocks[i + 1][3] >= self.min_supporting_variants and
                haplotype_blocks[i][2] != haplotype_blocks[i + 1][2]):
                
                is_recombinant = True
                
                end_pos1 = haplotype_blocks[i][1]
                start_pos2 = haplotype_blocks[i + 1][0]
                
                breakpoints.append(f"{end_pos1}-{start_pos2}")
        
        return is_recombinant, breakpoints

def parse_vcf(vcf_file):
    variants = {}
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
    parser.add_argument("-m", "--min_supporting_variants", type=int, default=2, help="Minimum number of supporting variants to consider a phase switch.")
    parser.add_argument("-w", "--workers", type=int, default=4, help="Number of worker processes to use.")
    args = parser.parse_args()

    variants = parse_vcf(args.vcf)
    results = process_bam(args.bam, variants, args.min_supporting_variants, num_workers=args.workers)

    for result in results:
        print(result)

if __name__ == "__main__":
    main()