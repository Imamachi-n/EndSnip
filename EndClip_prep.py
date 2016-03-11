#!/usr/bin/env python3

import sys
import os

###Custom modules
from module.utils_others import Input_config_file, now_time
from module.utils_extract_UTR_bed_file import Extract_3UTR_from_bed, Extract_gene_symbol_map_kfXref_file
from module.utils_extract_pA_sites import Merge_pA_site_infor

def main(config_file):
    #Greeting
    now_time('Beginning EndClip Prep run (v0.1.0)')
    print("-"*50)

    #Parse configure file
    now_time("Parsing configure file...")
    config_dict = Input_config_file(config_file)

    #Check configure file
    if not 'gene_gtf_file' in config_dict:
        sys.exit("ERROR: Gene GTF file does not exist...")
    if not 'gene_bed_file' in config_dict:
        sys.exit("ERROR: Gene BED file does not exist...")
    if not 'pA_site_file' in config_dict:
        sys.exit("ERROR: pA site file does not exist...")
    if not 'output_kfXref_file' in config_dict:
        sys.exit("ERROR: output_kfXref_file was not designated...")
    if not 'output_utr_file' in config_dict:
        sys.exit("ERROR: output_utr_file does not exist...")
    if not 'output_initial_UTR_database_file' in config_dict:
        sys.exit("ERROR: output_initial_UTR_database_file was not designated...")

    #Extract information
    gene_gtf_file = config_dict['gene_gtf_file']
    gene_bed_file = config_dict['gene_bed_file']
    pA_site_file = config_dict['pA_site_file']
    output_kfXref_file = config_dict['output_kfXref_file']
    output_utr_file = config_dict['output_utr_file']
    output_initial_UTR_database_file = config_dict['output_initial_UTR_database_file']

    #Prepare symbol_refid_map_file(kfXref file)
    '''
    now_time("Prepare kfXref file...")
    symbol_refid_map = Extract_gene_symbol_map_kfXref_file(gene_gtf_file)
    map_file = open(output_kfXref_file, 'w')
    for refid in symbol_refid_map.keys():
        symbol = symbol_refid_map[refid]
        print(refid, symbol, sep="\t", end="\n", file=map_file)
    '''

    #Make initial 3'UTR database
    now_time("Prepare initial 3'UTR database...")
    test = Extract_3UTR_from_bed(gene_bed_file, output_kfXref_file, output_utr_file)


    ##Prepare 3UTR database with pA site information
    #now_time("Prepare 3'UTR database with pA site information...")
    #temp_file = output_initial_UTR_database_file + '.tmp'
    #cmd = 'bedtools intersect -a %s -b %s -wa -wb > %s' % (output_utr_file, pA_site_file, temp_file)
    #os.system(cmd)

    #Merge_pA_site_infor(temp_file, output_initial_UTR_database_file)

    now_time("Completely finished!!")


if __name__ == '__main__':
    '''The identification of alternative polyadenylation sites 
       with RNA-seq, 3'-seq and poly(A) site database.
    '''
    main(sys.argv[1])
