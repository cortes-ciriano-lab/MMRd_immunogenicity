import math
import pysam
import timeit
import multiprocessing as mp
import argparse
import pandas as pd
import glob
import os
# import sys

def collect_result(result):
    if (result[1] != ''):
        VARIANTS = result[1]
        OUT = result[0]
        out = open(OUT,'w')
        out.write(VARIANTS)
        out.close()
        out_temp = result[0]


def concatenate_sort_temp_files_and_write(out_file, tmp_dir):
    # Get the file paths
    all_files = glob.glob(tmp_dir + '/*.BaseCounts.temp')
    
    # Load as panda files
    if (len(all_files) > 0):
        list_df = []
        
        for filename in all_files:
            df = pd.read_csv(filename, sep ='\t', index_col=None, header=None)
            list_df.append(df)
            os.remove(filename)
        
        # Concatenate and add header to df
        Merge = pd.concat(list_df, axis=0, ignore_index=True)

        # Add header
        Header = ['CHROM', 'POS', 'REF', 'ALT', 'DP', 'AF', 'A', 'C', 'T' , 'G', 'I', 'D', 'N', 'Other']
        Merge.columns = Header
     
        # Sort by chromosome and position
        final_sort = Merge.sort_values(['CHROM', 'POS'])
        
        # Write final file
        final_sort.to_csv (out_file, index = False, header=True, sep ='\t')
    else:
        # If temp files not found, print message
        print ('No temporary files found')

def BaseCount(LIST, REF_BASE):
    Bases=['A','C','T','G','N']
    
    # Dictinary with base counts
    NUCLEOTIDES = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, '+': 0, '-' : 0, 'N': 0, 'R' : 0, 'Other' : 0}
    
    # Update dictionary counts
    for x in LIST:
        if x.upper() in Bases:
            NUCLEOTIDES[x.upper()] += 1
        elif '-' in x:
            NUCLEOTIDES['-'] += 1
        elif '+' in x:
            NUCLEOTIDES['+'] += 1
        else:
            NUCLEOTIDES['Other'] += 1
    
    # Calculate Alternative count of the most frequent alternative allele
    Alternative = ['A','C','T','G','N','+', '-']
    Alternative = [x for x in Alternative if x != REF_BASE.upper()]
    
    alleles_by_freq = sorted(Alternative, key=lambda k: NUCLEOTIDES[k], reverse=True)
    ALT_ALLELE = alleles_by_freq[0]
    MAX_AC = NUCLEOTIDES[ALT_ALLELE]
    
    if (MAX_AC == 0):
        ALT_ALLELE = '.'

        
    return ([NUCLEOTIDES,ALT_ALLELE,MAX_AC])

#@profile
def run_interval(KEY, DICT, BAM, FASTA, MIN_COV, tmp_dir, BQ, MQ):
    
    interval = DICT[KEY]
    
    # Coordinates to analyse
    CHROM  =  interval[0]
    SITES = interval[1]
    
    # Get pileup read counts from coordinates
    bam = pysam.AlignmentFile(BAM)
    i = bam.pileup(CHROM, min(SITES)-1, max(SITES)+1, min_base_quality = BQ, min_mapping_quality = MQ, max_depth = 1000000, ignore_overlaps = False)
    
    # Load reference file. Mandatory to be done inside function to avoid overlap problems during multiprocessing
    inFasta = pysam.FastaFile(FASTA)

    # Run it for each position in pileup
    POSITIONS = []
    for p in i:
        POS=p.pos
        if POS in SITES:
            # Get reference base from fasta file
            ref_base = inFasta.fetch(CHROM, POS, POS+1)
            
            # Get coverage
            DP = p.get_num_aligned()
            
            # Run only if coverage is more than minimum (arg parameter)
            if (DP > MIN_COV and ref_base.upper() != 'N'):
                # Get pileup list
                PILEUP_LIST = p.get_query_sequences(mark_matches=True, add_indels=True)
                
                # Obtain base counts
                BASE_COUNTS, ALT_ALLELE, MAX_AC = BaseCount(PILEUP_LIST, ref_base)
                AF = float(MAX_AC) / DP

                
                # Print line only if minimum allele fraction and minimum alternative count is achieved
                POS_print = POS + 1
                    
                # Lines to print
                BC_A = str(BASE_COUNTS['A'])
                BC_C = str(BASE_COUNTS['C'])
                BC_T = str(BASE_COUNTS['T'])
                BC_G = str(BASE_COUNTS['G'])
                BC_I = str(BASE_COUNTS['+'])
                BC_D = str(BASE_COUNTS['-'])
                BC_N = str(BASE_COUNTS['N'])
                BC_Other = str(BASE_COUNTS['Other'])
                AF = str(round(AF, 4))
                
                LINE = [str(CHROM), str(POS_print), str(ref_base), str(ALT_ALLELE), str(DP), AF, BC_A, BC_C, BC_T, BC_G, BC_I, BC_D, BC_N, BC_Other]
                LINE = '\t'.join(LINE)
               
                # Save results in list of positions
                POSITIONS.append(LINE)
            
    inFasta.close()
    bam.close()
    
    # Return list of positions
    ID = '_'.join([str(CHROM), KEY])
    out_temp = tmp_dir + '/' + ID + '.BaseCounts.temp'
    return([out_temp,'\n'.join(POSITIONS)])

def split_genome(candidates, w_size):
    LIST_OF_SITES = {}
    with open(candidates, 'r') as cand:
        for line in cand:
            if not line.startswith('#'):
                line  =  line.rstrip('\n')
                LINE = line.split('\t')
                
                # Coordinates. They must be 0-based for pysam
                CHROM = str(LINE[0])
                START = int(LINE[1])
                SITE = int(LINE[2]) - 1 # In pysam, we used 0-based coordinates
                
                GENOME_WINDOW = math.floor(SITE / float(w_size))
                                
                # Key for dictionary for the pieces of the genome
                KEY = CHROM + '_' + str(GENOME_WINDOW)
                
                if KEY not in LIST_OF_SITES.keys():
                    LIST_OF_SITES[KEY] = {}
                    LIST_OF_SITES[KEY] = [CHROM, []]
                        
                LIST_OF_SITES[KEY][1].append(SITE)
    
    return (LIST_OF_SITES)

def initialize_parser():
    parser = argparse.ArgumentParser(description='Script to obtain base counts from bam based on a bed file')
    parser.add_argument('--bam', type=str, default=1, help='Tumor bam file to be analysed', required = True)
    parser.add_argument('--ref', type=str, default=1, help='Reference genome. *fai must be available in the same folder as reference', required = True)
    parser.add_argument('--out_file', help='Out file with potential somatic calls', required = True)
    parser.add_argument('--nprocs', type=int, default=1, help='Number of processes [Default = 1]',required=False)
    parser.add_argument('--bin', type=int, default=50000, help='Bin size for running the analysis', required = False)
    parser.add_argument('--bed', type=str, default='', help='Regions to focus the analysis. Three-column bed file', required = False)
    parser.add_argument('--min_cov', type=int, default = 1, help='Minimum coverage to consider the genomic site. [Default = 1]', required = False)
    parser.add_argument('--min_bq', type=int, default = 20, help='Minimum base quality permited for the base counts. [Default = 20]', required = False)
    parser.add_argument('--min_mq', type=int, default = 255, help='Minimum mapping quality required to analyse read. [Default = 255]', required = False)
    parser.add_argument('--tmp_dir', type=str, default = '.', help='Temporary folder for tmp files', required = False)

    return (parser)

def main():


    # 1. Arguments
    parser = initialize_parser()
    args = parser.parse_args()

    BAM = args.bam
    FASTA = args.ref
    CORE = args.nprocs
    out_file = args.out_file
    BIN = args.bin
    bed = args.bed
    MIN_COV = args.min_cov
    MIN_BQ = args.min_bq
    MIN_MQ = args.min_mq
    tmp_dir = args.tmp_dir


    # 2. Create bed file and windows
    WINDOW = split_genome(bed, BIN)
    
    # 3. Code to run in parallel all bins
    if (CORE > 1):
        pool = mp.Pool(CORE)
        
        # Step 3.1: Use loop to parallelize
        for KEY in WINDOW.keys():
            # This funtion writes in temp files the results
            pool.apply_async(run_interval, args=(KEY, WINDOW, BAM, FASTA, MIN_COV, tmp_dir, MIN_BQ, MIN_MQ), callback=collect_result)
                   
        # Step 3.2: Close Pool and let all the processes complete    
        pool.close()
        pool.join()
    else:
        for KEY in WINDOW:
            # This funtion writes in temp files the results
            collect_result(run_interval(KEY, WINDOW, BAM, FASTA, MIN_COV, tmp_dir, MIN_BQ, MIN_MQ))
            
    # 4. Write final file
    concatenate_sort_temp_files_and_write(out_file, tmp_dir)


if __name__ == '__main__':
    start = timeit.default_timer()
    main()
    print ('TIME \n')
    stop = timeit.default_timer()
    print (stop - start)
    


