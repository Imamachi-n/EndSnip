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

###Modules for testing
import matplotlib.pyplot as plt


def Load_Target_Wig_files(All_Wig_files, UTR_Annotation_file):
    UTR_events_dict = {}
    All_Samples_Total_depth = []

    #Load UTR Annotation file
    for line in open(UTR_Annotation_file, 'r'):
        fields = line.rstrip().split("\t")
        curr_chr = fields[0]
        region_start = fields[1]
        region_end = fields[2]
        name = fields[3]
        curr_strand = fields[5]
        UTR_pos = "%s:%s-%s" % (curr_chr, region_start, region_end)
        
        #Define 3'UTR Annotation regions
        #end_shift = int(round(abs(int(region_start) - int(region_end)) * 0.2)) # TODO: 2割の領域でいいか確認する。
        end_shift = 0
        if curr_strand == '+':
            region_end = str(int(region_end) - end_shift)
        elif curr_strand == '-':
            region_start = str(int(region_start) + end_shift)
        else:
            sys.exit("ERROR: Strand column in your UTR annotation file is wrong...")
        
        region_start = int(region_start) + 1
        region_end = int(region_end) - 1

        if (region_end - region_start) > 50: #Min 3UTR length(Default: 50bp)
            #UTR_events_dict => [chrom, start, end, strand, UTR_position(chrom:start-end)]
            UTR_events_dict[name] = [curr_chr, region_start, region_end, curr_strand, UTR_pos]
            # TODO: Isoformごとに判断する場合を考慮に入れる
            # 終止コドンが異なるケースでは、Isoformごとに判断する。
            # 3'UTR中でスプライシングを受けている場合も考慮する。

    #Load coverage for all samples
    All_samples_extracted_3UTR_coverage_dict = {}
    for curr_wig_file in All_Wig_files:
        curr_sample_All_chroms_coverage_dict = {}
        num_line = 0
        curr_sample_total_depth = 0
        for line in open(curr_wig_file, 'r'):
            if '#' in line and line[0:3] != 'chr':
                continue
            fields = line.strip().split("\t")
            chrom_name = fields[0]
            region_start = int(fields[1])
            region_end = int(fields[2])
            read_depth = int(float(fields[-1]))
            
            #Initialize coverage data in each chromosome
            if chrom_name not in curr_sample_All_chroms_coverage_dict:
                curr_sample_All_chroms_coverage_dict[chrom_name] = [[0],[0]] #[[region_site], [depth]]

            #Add coverage data in each region on each chromosome
            if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]: #if gap region exists
                curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start) #Region end => Region start
                curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)            #Read depth => 0
            curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_end)       #Region end
            curr_sample_All_chroms_coverage_dict[chrom_name][1].append(read_depth)       #Read depth
             
            #Total coverage and read count in each sample
            curr_sample_total_depth += read_depth * (region_end - region_start)
            num_line += 1

        #Collect total depth for each sample
        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0) #slicing your list(extracted_coverage): 185 line
        All_Samples_Total_depth.append(curr_sample_total_depth) #Collection of total depth for each sample
        
        #print(curr_sample_All_chroms_coverage_dict[chrom_name][0][1:20])
        #print(curr_sample_All_chroms_coverage_dict[chrom_name][1][1:20])
        #print(len(curr_sample_All_chroms_coverage_dict[chrom_name][0]))
        #print(len(curr_sample_All_chroms_coverage_dict[chrom_name][1]))

        now_time(curr_wig_file + ": Total depth were loaded.")

        #Define each depth for each 3'UTR
        for curr_3UTR_event_id in UTR_events_dict.keys():
            #Each transcript information
            curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
            curr_chr = curr_3UTR_structure[0] #Chromosome number

            #Call current chromosome from dictionary
            if curr_chr in curr_sample_All_chroms_coverage_dict.keys():
                #Region and Depth for current chromosome
                curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr]

                region_start = curr_3UTR_structure[1] #3'UTR region start
                region_end = curr_3UTR_structure[2] #3'UTR region end

                left_region_index = bisect(curr_chr_coverage[0], region_start) #Insertion site(index) of region start
                right_region_index = bisect(curr_chr_coverage[0], region_end) #Insertion site(index) of region end

                #List of depth(coverage) in 3'UTR region
                extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]

                #List of 3UTR region
                extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                extracted_3UTR_region.insert(0, region_start)
                extracted_3UTR_region.append(region_end)

                ###Example:
                ###chrom_site = [0,10,20,30,40,50,60,70,80,90,100]
                ###trx_exp =[0,5,5,5,5,5,1,1,0,0,0]
                ###bisect(chrom_site, 15)
                ###[2] => [0,10 | 20,30,40,50,60,70,80,90,100]
                ###bisect(chrom_site, 85)
                ###[9] => [0,10,20,30,40,50,60,70,80 | 90,100]
                ###chrom_site[2:9] => [20,30,40,50,60,70,80] (7 items) => [15,20,30,40,50,60,70,80,85] (8 items)
                ###trx_exp[2:9+1] => [5,5,5,5,1,1,0,0] (8 items)

                #Initiate current 3UTR event id in All_samples_extracted_3UTR_coverage_dict
                if not curr_3UTR_event_id in All_samples_extracted_3UTR_coverage_dict:
                    All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id] = []

                #Reserve Depth(Coverage) and 3UTR region information for each sample in dictonary
                All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id].append([extracted_coverage, extracted_3UTR_region]) 
                #Gene information(curr_3UTR_event_id) =>
                #[ [[Coverage list 1], [3UTR region list 1]], [[Coverage list 2], [3UTR region list 2]], ... , [[Coverage list N], [3UTR region list N]] ]

        now_time(curr_wig_file + ": Each depth for each 3'UTR was loaded.")

    #Reserve Depth(Coverage) and 3UTR region information for each 3UTR region | Total depth in samples | 3UTR region information
    return All_samples_extracted_3UTR_coverage_dict, np.array(All_Samples_Total_depth), UTR_events_dict

def Convert_wig_into_bp_coverage(extracted_coverage, extracted_3UTR_region, strand_info):
    #Initiate a list of 1bp-resolution coverage/chromosome site in 3UTR region for each gene
    UTR_length = extracted_3UTR_region[-1] - extracted_3UTR_region[0]
    bp_coverage = np.zeros(UTR_length)
    bp_chrom_site = list(range(extracted_3UTR_region[0],extracted_3UTR_region[-1])) #TODO:
    #print(len(bp_coverage))
    #print(len(bp_chrom_site))
    print(bp_chrom_site[0],bp_chrom_site[-1])

    #Start site in 3UTR region on chromosome
    relative_start = extracted_3UTR_region[0]

    #Make a list of 1bp-resolution coverage in 3UTR region from 5'end of last exon (0-base)
    for i in range(len(extracted_coverage)):
        curr_region_start = extracted_3UTR_region[i] - relative_start #start site in its transcript
        curr_region_end = extracted_3UTR_region[i+1] - relative_start #end site in its transcript
        bp_coverage[curr_region_start:curr_region_end] = extracted_coverage[i] #List of coverage

    #If stand is minus '-', reverse 1bp-resolution coverage list
    if strand_info == '-':
        bp_coverage = bp_coverage[::-1]

    #print(len(bp_coverage))
    #print(len(bp_chrom_site))
    bp_coverage = list(bp_coverage)
    return [bp_coverage,bp_chrom_site]

def De_Novo_3UTR_all_samples_bp_extimation(All_Samples_curr_3UTR_coverages, UTR_start, UTR_end, curr_strand, weight_for_second_coverage, Coverage_pPAS_cutoff):
    ###For each gene###
    #Parameter setting
    coverage_threshold = Coverage_pPAS_cutoff #Depth(Coverage) threshold in 5'end of last exon
    search_point_start = 200
    search_point_end = 200 #int(abs((UTR_end - UTR_start)) * 0.1) 
    # TODO: #3UTR長の「1割」でよいかどうか検討(3'-seqデータを元に3'endにpolyA siteが集中するケースが無いのかどうか検討する。)
    coverage_test_region = 100 #testing 0-100bp in 5'end of last exon
    least_3UTR_length = 500 #>=150bp 3UTR length are needed
    # TODO: least_3UTR_lengthはsearch_point_startとsearch_point_endの合計値よりも大きくなくてはいけない。あとでCriteriaを決める。

    #The number of samples
    num_samples = len(All_Samples_curr_3UTR_coverages)

    #Read coverage for each sample(List)
    Region_Coverages = [] #1bp coverage
    Region_mean_Coverages = [] #Mean of coverage
    Region_first_100_coverage_all_samples = [] #Mean of coverage in 100bp region (5'end last exon)

    #Prepare coverage information for each sample
    for i in range(num_samples):
        curr_Region_Coverage_raw = All_Samples_curr_3UTR_coverages[i] #Strand is reversed in load
        curr_Region_Coverage = curr_Region_Coverage_raw / weight_for_second_coverage[i] #bp_Coverage in 3UTR region for each sample

        Region_mean_Coverages.append(np.mean(curr_Region_Coverage_raw)) #Mean of coverage
        Region_Coverages.append(curr_Region_Coverage) #1bp coverage

        curr_first_100_coverage = np.mean(curr_Region_Coverage_raw[0:coverage_test_region]) #Mean of coverage in 100bp region (5'end last exon)
        Region_first_100_coverage_all_samples.append(curr_first_100_coverage) #List of mean coverage in 100bp region

    #Filtering coverage threshold(Default: >=5 coverage / >=150bp 3UTR length)
    if sum(np.array(Region_first_100_coverage_all_samples) >= coverage_threshold) >= num_samples and (UTR_end - UTR_start) >= least_3UTR_length:
        #Select 3UTR search region and make it's list
        if curr_strand == '+': #strand: +
            search_region = range(UTR_start + search_point_start, UTR_end - search_point_end + 1)
        elif curr_strand == '-': #strand: -
            search_region = range(UTR_end - search_point_start, UTR_start + search_point_end - 1, -1)

        #Define search region on 3UTR region(start/end)
        search_region_start = search_point_start #Default: 200bp
        search_region_end = UTR_end - UTR_start - search_point_end
        #print(search_region[0],search_region[-1])
        #print(UTR_start)
        #print(UTR_end)
        #print(search_region_start)
        #print(search_region_end)

        Mean_squared_error_list = []
        Estimated_3UTR_abundance_list = []

        #for each base-pair(bp)
        #print(list(range(search_region_start, search_region_end + 1)))
        for curr_point in range(search_region_start, search_region_end + 1):
            curr_search_point = curr_point #current base-pair on 3UTR region
            
            #Initiate a list of results for each sample
            All_samples_result = [ [], [], [] ]

            #testing coverage variance for each base-pair(bp) for each sample
            for curr_sample_region_coverage in Region_Coverages:
                Mean_Squared_error, Long_UTR_abun, Short_UTR_abun = Estimation_abundance(curr_sample_region_coverage, curr_search_point)
                #Output result
                All_samples_result[0].append(Mean_Squared_error)
                All_samples_result[1].append(Long_UTR_abun)
                All_samples_result[2].append(Short_UTR_abun)

            #DaPars Result (SampleA + SampleB) #TODO: ２種類のサンプルすべてのデータの平均値を出しているため、検出感度が落ちる原因になっている？
            Mean_Squared_error = np.mean(np.array(All_samples_result[0])) #Mean of coverage variance among samples
            Mean_squared_error_list.append(Mean_Squared_error) #Mean of coverage variance among each 3UTR set(short/long)
            Estimated_3UTR_abundance_list.append([All_samples_result[1], All_samples_result[2]]) #Long/short UTR abundance

        #TODO: Cross Point:CTRLとCFIm25ノックダウンでの各bpのcoverageの差分を求める。
        coverage_variance_CP = np.array(Region_Coverages[1] - Region_Coverages[0])
        coverage_variance_CP2 = coverage_variance_CP**2
        plt.plot(coverage_variance_CP2)
        filename = "data/output_CP_variance_" + test_name + ".png"
        plt.savefig(filename)
        plt.close()

        test_length = len(coverage_variance_CP)
        plt.plot(coverage_variance_CP)
        plt.plot(np.repeat([0],test_length))
        filename = "data/output_CP2_variance_" + test_name + ".png"
        plt.savefig(filename)
        plt.close()

        #Estimate min mean squared error
        if len(Mean_squared_error_list):
            #TODO: TEST: Mean_squared_error in 3'UTR region for PTEN, ELAVL1
            plt.plot(Mean_squared_error_list)
            #plt.show()
            filename = "data/output_variance_" + test_name + ".png"
            plt.savefig(filename)
            plt.close()

            #Identify index of min mean squared error in Mean_squared_error_list
            min_ele_index = Mean_squared_error_list.index(min(Mean_squared_error_list))

            #Result
            select_mean_squared_error = Mean_squared_error_list[min_ele_index] #Min mean squared error
            UTR_abundances = Estimated_3UTR_abundance_list[min_ele_index] #Long/short UTR abundance [Long UTR abundance, Short UTR abundance]
            selected_break_point = search_region[min_ele_index] #Chromosome site of break point

        else:
            select_mean_squared_error = 'NA'
            UTR_abundances = 'NA'
            selected_break_point = 'NA'
            print("WARNING: mean squared error list is empty...")

    else:
        select_mean_squared_error = 'NA'
        UTR_abundances = 'NA'
        selected_break_point = 'NA'

    return select_mean_squared_error, UTR_abundances, selected_break_point

def Estimation_abundance(Region_Coverage, break_point):
    Region_Coverage_long = Region_Coverage[break_point:] #List of coverage in long UTR region
    Region_Coverage_short = Region_Coverage[0:break_point] #List of coverage in short UTR region
    #print(Region_Coverage_short)
    #print(Region_Coverage_long)

    #Calculate mean of coverage(Long/Short UTR)
    Long_UTR_abun = np.mean(Region_Coverage_long)
    Short_UTR_abun = np.mean(Region_Coverage_short - Long_UTR_abun)
    if Short_UTR_abun < 0:
        Short_UTR_abun = 0
    #print(Long_UTR_abun)
    #print(Short_UTR_abun)

    #Calculate coverage variance
    Long_UTR_coverage_residual = Region_Coverage_short - Long_UTR_abun - Short_UTR_abun #List of coverage residual for short UTR region
    Short_UTR_coverage_residual = Region_Coverage_long - Long_UTR_abun #List of coverage residual for long UTR region
    Coverage_residual = Long_UTR_coverage_residual
    #print(Coverage_residual)
    Coverage_residual = np.append(Coverage_residual, Short_UTR_coverage_residual)
    
    Mean_Squared_error = np.mean(Coverage_residual**2) #Coverage variance for UTR region
    
    return Mean_Squared_error, Long_UTR_abun, Short_UTR_abun

def main():
    '''Input_configure_file
    Annotated_3UTR=hg19_refseq_extracted_3UTR.bed
    PolyA_site_infor=polyA_DB_hg19.bed
    Group1_Tophat_aligned_Wig=Condition_A_chrX.wig
    Group2_Tophat_aligned_Wig=Condition_B_chrX.wig
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
        sys.exit("ERROR: Group1_Tophat_aligned_Wig does not exist...")
    if not 'Group2_Tophat_aligned_Wig' in config_dict:
        sys.exit("ERROR: Group2_Tophat_aligned_Wig does not exist...")
    if not 'Output_directory' in config_dict:
        sys.exit("ERROR: Output_directory does not exist...")
    if not 'Annotated_3UTR' in config_dict:
        sys.exit("ERROR: Annotated_3UTR was not designated...")
    if not 'Output_result_file' in config_dict:
        sys.exit("ERROR: Output_result_file does not exist...")

    #File/Directory
    Group1_Tophat_aligned_file = config_dict['Group1_Tophat_aligned_Wig']
    Group2_Tophat_aligned_file = config_dict['Group2_Tophat_aligned_Wig']
    output_directory = config_dict['Output_directory']
    if output_directory[-1] != '/':
        output_directory += '/'
    Annotated_3UTR_file = config_dict['Annotated_3UTR']
    Output_result_file = config_dict['Output_result_file']

    #Parameter
    #TODO: リストを使って、ディクショナリのキーに含まれているかどうか判断。リスト内包表記などを利用。
    Num_least_in_group1_local = config_dict['Num_least_in_group1']
    Num_least_in_group2_local = config_dict['Num_least_in_group2']
    Coverage_cutoff_local = config_dict['Coverage_cutoff']
    FDR_cutoff_local = config_dict['FDR_cutoff']
    Fold_change_cutoff_local = config_dict['PDUI_cutoff']
    PDUI_cutoff_local = config_dict['Fold_change_cutoff']
    Coverage_pPAS_cutoff_local = config_dict['Coverage_pPAS_cutoff']

    #Group1_Tophat_aligned_file, Group2_Tophat_aligned_file, output_directory, Annotated_3UTR_file, PolyA_site_infor, Output_result_file, Num_least_in_group1_local, Num_least_in_group2_local, Coverage_cutoff_local, FDR_cutoff_local, Fold_change_cutoff_local, PDUI_cutoff_local, Coverage_pPAS_cutoff_local = parse_cfgfile(cfg_file)

    #Collect sample files
    num_group_1 = len(Group1_Tophat_aligned_file)
    num_group_2 = len(Group2_Tophat_aligned_file)

    All_Sample_files = Group1_Tophat_aligned_file[:]
    All_Sample_files.extend(Group2_Tophat_aligned_file)

    num_samples = len(All_Sample_files)

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
    Coverage_pPAS_cutoff = 5

    #Check parameters
    if Num_least_in_group1_local != '':
        Num_least_in_group1 = float(Num_least_in_group1_local)
    if Num_least_in_group2_local != '':
        Num_least_in_group2 = float(Num_least_in_group2_local)
    if Coverage_cutoff_local != '':
        Coverage_cutoff = float(Coverage_cutoff_local)
    if FDR_cutoff_local != '':
        FDR_cutoff = float(FDR_cutoff_local)
    if Fold_change_cutoff_local != '':
        Fold_change_cutoff = float(Fold_change_cutoff_local)
    if PDUI_cutoff_local != '':
        PDUI_cutoff = float(PDUI_cutoff_local)
    if Coverage_pPAS_cutoff_local != '':
        Coverage_pPAS_cutoff = float(Coverage_pPAS_cutoff_local)

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
    first_line = ['Gene','fit_value','Predicted_Proximal_APA','loci']
    for i in range(num_group_1):
        curr_long_exp = 'A_%s_long_exp' % str(i+1)
        curr_short_exp = 'A_%s_short_exp' % str(i+1)
        curr_ratio = 'A_%s_PDUI' % str(i+1)
        first_line.extend([curr_long_exp, curr_short_exp, curr_ratio])
    for i in range(num_group_2):
        curr_long_exp = 'A_%s_long_exp' % str(i+1)
        curr_short_exp = 'A_%s_short_exp' % str(i+1)
        curr_ratio = 'A_%s_PDUI' % str(i+1)
        first_line.extend([curr_long_exp, curr_short_exp, curr_ratio])
    first_line.append('PDUI_Group_diff')

    print("\t".join(first_line), end="\n", file=Output_result)

    #Test APA event for each 3UTR
    now_time("Testing APA events for each 3UTR region...")
    for curr_3UTR_id in UTR_events_dict:
        #3UTR region information for each gene
        curr_3UTR_structure = UTR_events_dict[curr_3UTR_id]
        region_start = curr_3UTR_structure[1] #region start
        region_end = curr_3UTR_structure[2] #region end
        curr_strand = curr_3UTR_structure[-2] #strand
        UTR_pos = curr_3UTR_structure[-1] #UTR position information

        #If gene names exist in coverage dict(for each gene)
        if curr_3UTR_id in All_samples_Target_3UTR_coverages:
            #3UTR coverage for each gene
            curr_3UTR_coverage_wig = All_samples_Target_3UTR_coverages[curr_3UTR_id] #List of 3UTR coverage for each sample
            curr_3UTR_all_samples_bp_coverage = []
            curr_3UTR_all_samples_bp_chrom_site = []
            for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig: #3UTR coverage for each sample
                bp_resolution_data = Convert_wig_into_bp_coverage(curr_sample_curr_3UTR_coverage_wig[0], #List of coverage
                                                                  curr_sample_curr_3UTR_coverage_wig[1], #List of 3UTR region
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
            select_mean_squared_error, selected_break_point, UTR_abundance = De_Novo_3UTR_all_samples_bp_extimation(curr_3UTR_all_samples_bp_coverage,
                                                                                                                    region_start,
                                                                                                                    region_end,
                                                                                                                    curr_strand,
                                                                                                                    All_sample_coverage_weights,
                                                                                                                    Coverage_pPAS_cutoff) 
            

if __name__ == '__main__':
    '''The identification of alternative polyadenylation sites 
       with RNA-seq, 3'-seq and poly(A) site database.
    '''
    main()
