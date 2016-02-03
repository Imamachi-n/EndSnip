#!/usr/bin/env python3
"""
EndClip: Comprehensive analysis of Alternative PolyAdenylation Site.
Created by Naoto Imamachi on 2016-01-06.
Updated and maintained by Naoto Imamachi since Jan 2016.

This script was modified from Dapars algorithm (https://github.com/ZhengXia/DaPars).
License: GPL v2

Usage:
  EndClip_main.py <configure_file>

"""

import sys
import os
from datetime import datetime

from bisect import bisect
import numpy as np

###Custom modules
from module.utils_others import now_time, Input_config_file
from module.core_load_wig_files import Load_Target_Wig_files, Convert_wig_into_bp_coverage
from module.core_coverage_comparison import De_Novo_3UTR_all_samples_bp_extimation, Estimation_abundance, coverage_comparison_with_pA_site, de_novo_coverage_comparison_with_windows

###Modules for testing
import matplotlib.pyplot as plt

def main():
    '''Input_configure_file
    Annotated_3UTR=data/hg19_refseq_extracted_3UTR_PTEN_ELAVL1.bed
    PolyA_site_infor=data/polyA_DB_hg19.bed
    Group1_Tophat_aligned_Wig=data/siCTRL_S_accepted_hits_PTEN_ELAVL1.bam.wig
    Group2_Tophat_aligned_Wig=data/siCFIm25-1_accepted_hits_PTEN_ELAVL1.bam.wig
    Output_directory=DaPars_Test_data/
    Output_result_file=DaPars_Test_data

    #Parameters
    Num_least_in_group1=1
    Num_least_in_group2=1
    Coverage_cutoff=30
    FDR_cutoff=0.05
    PDUI_cutoff=0.5
    Fold_change_cutoff=0.59
    '''
    now_time("Beginning EndClip run (v0.1.0)")
    print("-"*50)
    if len(sys.argv) == 1:
        sys.exit("ERROR: Please provide the configure file...")
    cfg_file = sys.argv[1]

    now_time("Parsing configure file...")
    config_dict = Input_config_file(cfg_file)

    #Check configure file
    if not 'Group1_Tophat_aligned_Wig' in config_dict:
        sys.exit("ERROR: No Tophat aligned BAM file for group 1...")
    if not 'Group2_Tophat_aligned_Wig' in config_dict:
        sys.exit("ERROR: No Tophat aligned BAM file for group 2...")
    if not 'Output_directory' in config_dict:
        sys.exit("ERROR: No output directory...")
    if not 'Annotated_3UTR' in config_dict:
        sys.exit("ERROR: No annotated 3'UTR file...")
    if not 'Output_result_file' in config_dict:
        sys.exit("ERROR: No result file name...")

    #File/Directory
    Group1_Tophat_aligned_file = config_dict['Group1_Tophat_aligned_Wig'].split(',')
    Group2_Tophat_aligned_file = config_dict['Group2_Tophat_aligned_Wig'].split(',')
    output_directory = config_dict['Output_directory']
    if output_directory[-1] != '/':
        output_directory += '/'
    Annotated_3UTR_file = config_dict['Annotated_3UTR']
    Output_result_file = config_dict['Output_result_file']

    #Default parameters
    global Num_least_in_group1
    global Num_least_in_group2
    global Coverage_cutoff
    global FDR_cutoff
    global Fold_change_cutoff
    global PDUI_cutoff
    global Coverage_pPAS_cutoff

    Num_least_in_group1 = 1
    Num_least_in_group2 = 1
    Coverage_cutoff = 30
    FDR_cutoff = 0.05
    Fold_change_cutoff = 0.59 #1.5-fold change
    PDUI_cutoff = 0.2
    Coverage_pPAS_cutoff = 5.0
    
    #Check parameters
    if not 'Num_least_in_group1' in config_dict:
        print("  Num_least_in_group1: Default parameter(1) was designated.")
    else:
        Num_least_in_group1 = float(config_dict['Num_least_in_group1'])

    if not 'Num_least_in_group2' in config_dict:
        print("  Num_least_in_group2: Default parameter(1) was designated.")
    else:
        Num_least_in_group2 = float(config_dict['Num_least_in_group2'])

    if not 'Coverage_cutoff' in config_dict:
        print("  Coverage_cutoff: Default parameter(30) was designated.")
    else:
        Coverage_cutoff = float(config_dict['Coverage_cutoff'])

    if not 'FDR_cutoff' in config_dict:
        print("  FDR_cutoff: Default parameter(0.05) was designated.")
    else:
        FDR_cutoff = float(config_dict['FDR_cutoff'])
    
    if not 'Fold_change_cutoff' in config_dict:
        print("  Fold_change_cutoff: Default parameter(0.59[log2]/1.5-fold) was designated.")
    else:
        Fold_change_cutoff = config_dict['Fold_change_cutoff']

    if not 'PDUI_cutoff' in config_dict:
        print("  PDUI_cutoff: Default parameter(0.2) was designated.")
    else:
        PDUI_cutoff = float(config_dict['PDUI_cutoff'])

    if not 'Coverage_pPAS_cutoff' in config_dict:
        print("  Coverage_pPAS_cutoff: Default parameter(5.0) was designated.")
    else:
        Coverage_pPAS_cutoff = float(config_dict['Coverage_pPAS_cutoff'])

    #Collect sample files
    num_group_1 = len(Group1_Tophat_aligned_file)
    num_group_2 = len(Group2_Tophat_aligned_file)

    All_Sample_files = Group1_Tophat_aligned_file[:]
    All_Sample_files.extend(Group2_Tophat_aligned_file)

    num_samples = len(All_Sample_files)

    #Prepare output directory
    d = os.path.dirname(output_directory)
    if not os.path.exists(d):
        os.makedirs(d)
    
    #Prepare temp directory
    temp_dir = d + '/tmp/'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    Output_all_prediction_file = output_directory + Output_result_file + '_result_temp.txt'
    Output_result = open(Output_all_prediction_file, 'w')

    #Load coverage
    now_time("Loading coverage...")
    All_samples_Target_3UTR_coverages, All_samples_sequencing_depths, UTR_events_dict = Load_Target_Wig_files(All_Sample_files, Annotated_3UTR_file)

    #Depth(Coverage) weight for each sample
    All_sample_coverage_weights = All_samples_sequencing_depths / np.mean(All_samples_sequencing_depths)
    now_time("Loading coverage finished.")

    #Prepare header information for output file
    first_line = ['Gene','Predicted_Proximal_APA','loci']
    for i in range(num_group_1):
        curr_long_exp = 'A_%s_long_exp' % str(i+1)
        curr_short_exp = 'A_%s_short_exp' % str(i+1)
        curr_ratio = 'A_%s_PDUI' % str(i+1)
        first_line.extend([curr_long_exp, curr_short_exp, curr_ratio])
    for i in range(num_group_2):
        curr_long_exp = 'B_%s_long_exp' % str(i+1)
        curr_short_exp = 'B_%s_short_exp' % str(i+1)
        curr_ratio = 'B_%s_PDUI' % str(i+1)
        first_line.extend([curr_long_exp, curr_short_exp, curr_ratio])
    first_line.extend(['A_PDUI_mean','B_PDUI_mean','PDUI_Group_diff[B_PDUI_mean - A_PDUI_mean]','Fold-change[A_PDUI_mean/B_PDUI_mean]','Fold-change[B_short_exp_mean / A_short_exp_mean]'])

    print("\t".join(first_line), end="\n", file=Output_result)

    #Test APA event for each 3UTR
    now_time("Testing APA events for each 3UTR region...")
    test_wig_output_file1 = open('test_bedgraph_file_ELAVL1_PTEN_siCTRL.bg', 'w') ###TEST:
    test_wig_output_file2 = open('test_bedgraph_file_ELAVL1_PTEN_siCFIm25.bg', 'w') ###TEST:
    print('track type=bedGraph name=EndClip_test_ELAVL1_PTEN_siCTRL description=EndClip_test_ELAVL1_PTEN_siCTRL visibility=2 maxHeightPixels=40:40:20', end="\n", file=test_wig_output_file1)
    print('track type=bedGraph name=EndClip_test_ELAVL1_PTEN_siCFIm25 description=EndClip_test_ELAVL1_PTEN_siCFIm25 visibility=2 maxHeightPixels=40:40:20', end="\n", file=test_wig_output_file2)

    for curr_3UTR_id in UTR_events_dict:
        #3UTR region information for each gene
        curr_3UTR_structure = UTR_events_dict[curr_3UTR_id]
        chrom = curr_3UTR_structure[0]
        region_start = curr_3UTR_structure[1] #region start
        region_end = curr_3UTR_structure[2] #region end
        curr_strand = curr_3UTR_structure[3] #strand
        UTR_pos = curr_3UTR_structure[4] #UTR position information
        pA_site = curr_3UTR_structure[5].split('|') #pA_site list
        pA_site = list(map(int,pA_site))

        #If gene names exist in coverage dict(for each gene)
        if curr_3UTR_id in All_samples_Target_3UTR_coverages:
            #3UTR coverage for each gene
            curr_3UTR_coverage_wig = All_samples_Target_3UTR_coverages[curr_3UTR_id] #List of 3UTR coverage for each sample
            curr_3UTR_all_samples_bp_coverage = []
            curr_3UTR_all_samples_bp_chrom_site = []
            for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig: #3UTR coverage for each sample
                bp_resolution_data = Convert_wig_into_bp_coverage(curr_sample_curr_3UTR_coverage_wig[0], #List of coverage
                                                                  curr_sample_curr_3UTR_coverage_wig[1], #List of 3UTR region(1-base)
                                                                  curr_strand) #strand
                #test = Convert_wig_into_bp_coverage(curr_sample_curr_3UTR_coverage_wig[0],curr_sample_curr_3UTR_coverage_wig[1],curr_strand) #test
                curr_3UTR_curr_samples_bp_coverage = bp_resolution_data[0]
                curr_3UTR_curr_samples_bp_chrom_site = bp_resolution_data[1]
                curr_3UTR_all_samples_bp_coverage.append(curr_3UTR_curr_samples_bp_coverage) #List of bp_coverage for each sample
                curr_3UTR_all_samples_bp_chrom_site.append(curr_3UTR_curr_samples_bp_chrom_site) #List of bp chromosome site for each sample
                
                #TODO: TEST: Coverage in 3'UTR region for PTEN, ELAVL1
                plt.plot(curr_3UTR_curr_samples_bp_coverage)
                #plt.show()
                global test_name
                test_name = curr_3UTR_id.split('|')[1]
                filename = "data/output_coverage_" + test_name + ".png"
                plt.savefig(filename)

            #TODO: TEST: Coverage in 3'UTR region for PTEN, ELAVL1
            plt.close()
            
            #De novo identification of APA event for each 3UTR region
            curr_3UTR_all_samples_bp_coverage = np.array(curr_3UTR_all_samples_bp_coverage)
            #select_mean_squared_error, selected_break_point, UTR_abundance = De_Novo_3UTR_all_samples_bp_extimation(curr_3UTR_all_samples_bp_coverage,
            #                                                                                                        region_start,
            #                                                                                                        region_end,
            #                                                                                                        curr_strand,
            #                                                                                                        All_sample_coverage_weights,
            #                                                                                                        Coverage_pPAS_cutoff,
            #                                                                                                        test_name) 
            #coverage_comparison_with_pA_site(curr_3UTR_all_samples_bp_coverage, curr_3UTR_all_samples_bp_chrom_site, region_start, region_end, curr_strand, All_sample_coverage_weights, Coverage_pPAS_cutoff, pA_site,test_name)
            de_novo_coverage_comparison_with_windows(curr_3UTR_all_samples_bp_coverage, curr_3UTR_all_samples_bp_chrom_site, region_start, region_end, curr_strand, All_sample_coverage_weights, Coverage_pPAS_cutoff, pA_site,test_name, chrom, test_wig_output_file1, test_wig_output_file2, Output_result, num_group_1, num_group_2, curr_3UTR_id, UTR_pos)

    now_time("Completely finished!!")

if __name__ == '__main__':
    '''The identification of alternative polyadenylation sites 
       with RNA-seq, 3'-seq and poly(A) site database.
    '''
    main()
