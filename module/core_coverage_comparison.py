#!/usr/bin/env python3

import math
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

from module.utils_coverage_comparison import *

###MAIN###
def de_novo_coverage_comparison_with_windows(curr_3UTR_all_samples_bp_coverage, curr_3UTR_all_samples_bp_chrom_site, UTR_start, UTR_end, curr_strand, weight_for_second_coverage, Coverage_pPAS_cutoff, pA_site, test_name, chrom, test_wig_output_file1, test_wig_output_file2, Output_result, num_group_1, num_group_2, curr_3UTR_id, UTR_pos):
    #curr_3UTR_all_samples_bp_coverage: 
    #[ [[Coverage list 1], [3UTR region list 1]], [[Coverage list 2], [3UTR region list 2]], ... , [[Coverage list N], [3UTR region list N]] ]
    ###For each gene###
    #Parameter setting
    coverage_threshold = Coverage_pPAS_cutoff #Depth(Coverage) threshold in 5'end of last exon
    search_point_start = 100
    search_point_end = 100
    coverage_test_region = 100 #testing 0-100bp in 5'end of last exon
    least_3UTR_length = 500 #>=150bp 3'UTR length are needed
    least_search_region_coverage = Coverage_pPAS_cutoff #>=5 coverage are needed in search region
    search_region_distance = 200 #Sub-region length in 3'UTR region #TODO: どうして200bpである必要があるのか理由付けがまるでない。

    line_write = []

    #The number of samples
    num_samples = len(curr_3UTR_all_samples_bp_coverage)

    #Read coverage for each sample(List)
    Region_Coverages = [] #1bp coverage (Normalized)
    Region_mean_Coverages = [] #Mean of coverage (Raw)
    Region_first_100_coverage_all_samples = [] #Mean of coverage in 100bp region (5'end last exon) (Raw)

    #Prepare coverage information for each sample
    for i in range(num_samples):
        curr_Region_Coverage_raw = curr_3UTR_all_samples_bp_coverage[i] #Strand is reversed in load
        curr_Region_Coverage = curr_Region_Coverage_raw / weight_for_second_coverage[i] #bp_Coverage in 3UTR region for each sample

        Region_mean_Coverages.append(np.mean(curr_Region_Coverage_raw)) #Mean of coverage
        Region_Coverages.append(curr_Region_Coverage) #1bp coverage

        curr_first_100_coverage = np.mean(curr_Region_Coverage_raw[0:coverage_test_region]) #Mean of coverage in 100bp region (5'end last exon)
        Region_first_100_coverage_all_samples.append(curr_first_100_coverage) #List of mean coverage in 100bp region

    #Filtering coverage threshold(Default: >=5 coverage / >=150bp 3UTR length)
    if sum(np.array(Region_first_100_coverage_all_samples) >= coverage_threshold) >= num_samples and (UTR_end - UTR_start) >= least_3UTR_length:
        #Prepare UTR search region
        UTR_start_site = int(UTR_start) #UTR_start is '1-base'.
        UTR_end_site = int(UTR_end)
        UTR_end_list = pA_site

        UTR_search_region = De_novo_UTR_search_region(UTR_start, UTR_end, curr_strand, search_region_distance)
        print(UTR_search_region)

        #Initialize break point infor list
        gathered_break_point = []

        #Estimate_each_sample
        flg = 0 #TODO: TEST:
        for curr_3UTR_curr_sample_bp_coverage in curr_3UTR_all_samples_bp_coverage:
            ###STEP1:
            break_point_candidates = Estimate_break_point(curr_3UTR_curr_sample_bp_coverage, UTR_search_region, UTR_start_site, UTR_end_site, curr_strand, search_point_start, search_point_end, least_search_region_coverage, test_name, flg)
            print(break_point_candidates)

            ##STEP2:
            #ReEstimate break points
            updated_UTR_search_region = Define_UTR_search_region(UTR_start_site, UTR_end_site, curr_strand, break_point_candidates)
            print(updated_UTR_search_region)
            updated_break_point_candidates = Estimate_break_point(curr_3UTR_curr_sample_bp_coverage, updated_UTR_search_region, UTR_start_site, UTR_end_site, curr_strand, search_point_start, search_point_end, least_search_region_coverage, test_name, flg)
            print("Updated_break_point_candidates")
            print(updated_break_point_candidates)

            ###STEP3: 
            #Near <=200bp break points => ReEstimate
            updated_break_point_candidates_testing = []
            test_dataset = updated_break_point_candidates[:]
            if curr_strand == '+':
                test_dataset.append(UTR_end)
                updated_break_point_candidates_testing = [test_dataset[x+1]-test_dataset[x]+1 for x in range(len(test_dataset)-1)]
                print(updated_break_point_candidates_testing)
            elif curr_strand == '-':
                test_dataset.append(UTR_start)
                updated_break_point_candidates_testing = [abs(test_dataset[x+1]-test_dataset[x])+1 for x in range(len(test_dataset)-1)]
                print(updated_break_point_candidates_testing)
            
            #print(sum(np.array(updated_break_point_candidates_testing) >= 200))
            #print(len(updated_break_point_candidates))
            if sum(np.array(updated_break_point_candidates_testing) >= 200) < len(updated_break_point_candidates):
                print("Retry...")
                #if curr_strand == '-':
                #    updated_break_point_candidates = updated_break_point_candidates[::-1]
                updated_UTR_search_region = Define_UTR_search_region(UTR_start_site, UTR_end_site, curr_strand, updated_break_point_candidates)
                print(updated_UTR_search_region)
                updated_break_point_candidates = Estimate_break_point(curr_3UTR_curr_sample_bp_coverage, updated_UTR_search_region, UTR_start_site, UTR_end_site, curr_strand, search_point_start, search_point_end, least_search_region_coverage, test_name, flg)
                print("Retry: Updated_break_point_candidates")
                print(updated_break_point_candidates)
            else: #TEST:
                pass
                #print("If Retry...")
                #updated_UTR_search_region = Define_UTR_search_region(UTR_start_site, UTR_end_site, curr_strand, updated_break_point_candidates)
                #print(updated_UTR_search_region)
                #updated_break_point_candidates = Estimate_break_point(curr_3UTR_curr_sample_bp_coverage, updated_UTR_search_region, UTR_start_site, UTR_end_site, curr_strand, search_point_start, search_point_end, least_search_region_coverage, test_name, flg)
                #print("Retry??: Updated_break_point_candidates")
                #print(updated_break_point_candidates)

            #Save updated break points
            gathered_break_point.extend(updated_break_point_candidates)

            flg = 1 #TODO: TEST:
        
        #DIFF: Comarison of Control vs CFIm25 knockdown samples
        #Sort gathered_break_points
        gathered_break_point = sorted(gathered_break_point)
        if curr_strand == '-':
            gathered_break_point = gathered_break_point[::-1]
        print(gathered_break_point)

        #Prepare filtered UTR search regions
        gathered_UTR_search_region = Define_UTR_search_region(UTR_start_site, UTR_end_site, curr_strand, gathered_break_point)
        print('Gathered_UTR_search_region')
        print(gathered_UTR_search_region)

        #Estimate merged region coverage in all samples
        merged_region_coverages = []
        Region_Coverages_array = np.array(Region_Coverages)
        merged_region_coverages = np.sum(Region_Coverages_array, axis=0)
        print(merged_region_coverages)

        break_point_for_diff = Estimate_break_point(merged_region_coverages, gathered_UTR_search_region, UTR_start_site, UTR_end_site, curr_strand, search_point_start, search_point_end, least_search_region_coverage, test_name, flg)
        print(break_point_for_diff)
        if curr_strand == '+':
            line_write = [curr_3UTR_id, '|'.join(list(map(str,break_point_for_diff))),UTR_pos]
        elif curr_strand == '-':
            line_write = [curr_3UTR_id, '|'.join(list(map(str,break_point_for_diff[::-1]))),UTR_pos]

        #Prepare each UTR isoform length
        UTR_isoform_length = []
        if curr_strand == '+':
            break_point_for_diff_test = break_point_for_diff[:]
            break_point_for_diff_test.append(UTR_end)
            for end_point in break_point_for_diff_test:
                curr_UTR_length = end_point - UTR_start
                UTR_isoform_length.append(curr_UTR_length)
        elif curr_strand == '-':
            break_point_for_diff_test = break_point_for_diff[:]
            break_point_for_diff_test.append(UTR_start)
            for end_point in break_point_for_diff_test:
                curr_UTR_length = UTR_end - end_point
                UTR_isoform_length.append(curr_UTR_length)
        print(UTR_isoform_length)

        #Estimate multi-UTR coverage level for each sample
        if curr_strand == '+':
            break_point_for_diff.insert(0, UTR_start)
            break_point_for_diff.append(UTR_end)
        elif curr_strand == '-':
            break_point_for_diff.insert(0, UTR_end)
            break_point_for_diff.append(UTR_start)

        #Estimate the coverage of each UTR isoform for each sample
        flg_test = 0
        multi_UTR_coverage = []
        Each_UTR_coverage = []
        Each_UTR_coverage_percentage = []
        for curr_3UTR_curr_sample_bp_coverage in curr_3UTR_all_samples_bp_coverage:
            coverage_infor = Estimate_UTR_isoform_expression(curr_3UTR_curr_sample_bp_coverage, break_point_for_diff, curr_strand, chrom, test_wig_output_file1, test_wig_output_file2, flg_test)
            multi_UTR_coverage.append(coverage_infor[0]) #[[75.0, 39.0, 14.0, 2.0], [228.0, 5.0, 1.0, 1.0]]
            Each_UTR_coverage.append(coverage_infor[1]) #[[36.0, 25.0, 12.0, 2.0], [223.0, 4.0, 0.0, 1.0]]
            Each_UTR_coverage_percentage.append(coverage_infor[2]) #[[0.47999999999999998, 0.33333333333333331, 0.16, 0.026666666666666668], [0.97807017543859653, 0.017543859649122806, 0.0, 0.0043859649122807015]]
            flg_test = 1

        multi_UTR_coverage_sample1 = np.array(multi_UTR_coverage[:num_group_1])
        multi_UTR_coverage_sample1_sum = np.sum(multi_UTR_coverage_sample1, axis=0)
        multi_UTR_coverage_sample2 = multi_UTR_coverage[num_group_1:]
        multi_UTR_coverage_sample2_sum = np.sum(multi_UTR_coverage_sample2, axis=0)

        Each_UTR_coverage_sample1 = np.array(Each_UTR_coverage[:num_group_1])
        Each_UTR_coverage_sample1_sum = np.sum(Each_UTR_coverage_sample1, axis=0)
        Each_UTR_coverage_sample2 = np.array(Each_UTR_coverage[num_group_1:])
        Each_UTR_coverage_sample2_sum = np.sum(Each_UTR_coverage_sample2, axis=0)

        Each_UTR_coverage_percentage_sample1 = np.array(Each_UTR_coverage_percentage[:num_group_1])
        Each_UTR_coverage_percentage_sample1_sum = np.sum(Each_UTR_coverage_percentage_sample1, axis=0)
        Each_UTR_coverage_percentage_sample2 = np.array(Each_UTR_coverage_percentage[num_group_1:])
        Each_UTR_coverage_percentage_sample2_sum = np.sum(Each_UTR_coverage_percentage_sample2, axis=0)

        #Each_UTR_coverage_sub = np.array(Each_UTR_coverage[1]) - np.array(Each_UTR_coverage[0])
        #Each_UTR_coverage_percentage_sub = np.array(Each_UTR_coverage_percentage[1]) - np.array(Each_UTR_coverage_percentage[0])

        print(multi_UTR_coverage)
        print(multi_UTR_coverage_sample1_sum)
        print(multi_UTR_coverage_sample2_sum)
        print(Each_UTR_coverage)
        print(Each_UTR_coverage_sample1_sum)
        print(Each_UTR_coverage_sample2_sum)
        print(Each_UTR_coverage_percentage)
        print(Each_UTR_coverage_percentage_sample1_sum)
        print(Each_UTR_coverage_percentage_sample2_sum)
        #print(Each_UTR_coverage_sub)
        #print(Each_UTR_coverage_percentage_sub)

        added_line_write = Estimate_PDUI_score(Each_UTR_coverage, Each_UTR_coverage_percentage, UTR_isoform_length, num_group_1, num_group_2)
        line_write.extend(added_line_write)

        #Result reporting
        print("\t".join(line_write), end="\n", file=Output_result)


########################################################################
###OLD_SCRIPTS##########################################################
########################################################################

def De_Novo_3UTR_all_samples_bp_extimation(All_Samples_curr_3UTR_coverages, UTR_start, UTR_end, curr_strand, weight_for_second_coverage, Coverage_pPAS_cutoff, test_name):
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




def coverage_comparison_with_pA_site(curr_3UTR_all_samples_bp_coverage, curr_3UTR_all_samples_bp_chrom_site, UTR_start, UTR_end, curr_strand, weight_for_second_coverage, Coverage_pPAS_cutoff, pA_site, test_name):
    ###For each gene###
    #Parameter setting
    coverage_threshold = Coverage_pPAS_cutoff #Depth(Coverage) threshold in 5'end of last exon
    search_point_start = 100
    search_point_end = 100
    coverage_test_region = 100 #testing 0-100bp in 5'end of last exon
    least_3UTR_length = 500 #>=150bp 3'UTR length are needed
    least_search_region_coverage = Coverage_pPAS_cutoff #>=5 coverage are needed in search region

    #The number of samples
    num_samples = len(curr_3UTR_all_samples_bp_coverage)

    #Read coverage for each sample(List)
    Region_Coverages = [] #1bp coverage (Normalized)
    Region_mean_Coverages = [] #Mean of coverage (Raw)
    Region_first_100_coverage_all_samples = [] #Mean of coverage in 100bp region (5'end last exon) (Raw)

    #Prepare coverage information for each sample
    for i in range(num_samples):
        curr_Region_Coverage_raw = curr_3UTR_all_samples_bp_coverage[i] #Strand is reversed in load
        curr_Region_Coverage = curr_Region_Coverage_raw / weight_for_second_coverage[i] #bp_Coverage in 3UTR region for each sample

        Region_mean_Coverages.append(np.mean(curr_Region_Coverage_raw)) #Mean of coverage
        Region_Coverages.append(curr_Region_Coverage) #1bp coverage

        curr_first_100_coverage = np.mean(curr_Region_Coverage_raw[0:coverage_test_region]) #Mean of coverage in 100bp region (5'end last exon)
        Region_first_100_coverage_all_samples.append(curr_first_100_coverage) #List of mean coverage in 100bp region

    #Filtering coverage threshold(Default: >=5 coverage / >=150bp 3UTR length)
    if sum(np.array(Region_first_100_coverage_all_samples) >= coverage_threshold) >= num_samples and (UTR_end - UTR_start) >= least_3UTR_length:
        #Prepare UTR search region
        UTR_start_site = int(UTR_start) #UTR_start is '1-base'.
        UTR_end_site = int(UTR_end)
        UTR_end_list = pA_site

        #Define UTR search region
        UTR_search_region = Define_UTR_search_region(UTR_start_site, UTR_end_site, curr_strand, UTR_end_list)
        print(UTR_search_region)

        #Estimate_each_sample
        flg = 0 #TODO: TEST:
        for curr_3UTR_curr_sample_bp_coverage in curr_3UTR_all_samples_bp_coverage:
            Estimate_break_point(curr_3UTR_curr_sample_bp_coverage, UTR_search_region, UTR_start_site, UTR_end_site, curr_strand, search_point_start, search_point_end, least_search_region_coverage, test_name, flg)
            flg = 1 #TODO: TEST:


def UTR_end_list_filtering_forward(UTR_end_list):
    #Remove pA_site near the other pA_site(Default: >=150bp)
    ###TEST(Strand:-): UTR_end_list = [1,200,400,450,600,800,850,890,1000]
    ###*      *      -   *        *    -   -   -    *  check
    ###|      |      |   |        |    |   |   |    |
    ###1     200    400 450      600  800 850 890 1000 chrom_site

    ###TEST(Strand:+): UTR_end_list = [1000,800,610,550,500,300,1]
    ###  *       *      -   -   *      *        * check
    ###  |       |      |   |   |      |        |
    ###1000     800    610 550 500    300       1 chrom_site
    UTR_end_list_filtered = [UTR_end_list[x] for x in range(len(UTR_end_list)) if UTR_end_list[x] == UTR_end_list[-2] or UTR_end_list[x] == UTR_end_list[-1] or abs(UTR_end_list[x+1] - UTR_end_list[x]) + 1 >= 200] #TODO: 150-200bpのどの長さがよいか検討。
    # UTR_start(UTR_end_list[-1]) & first pA_site(UTR_end_list[-2]) are absolutely reserved in UTR_end_list_filtered list.

    #Reverse a list (if curr_stand is '-')
    if curr_strand == '+':
        UTR_end_list_filtered = UTR_end_list_filtered[::-1]
    elif curr_strand == '-':
        UTR_end_list_filtered = UTR_end_list_filtered[::-1]

    #Make triple dataset
    #UTR_end_list = [10,20,30,40,45]
    #UTR_start_site = 0
    #UTR_end_site = 50
    #UTR_search_region: [[0, 10, 20], [10, 20, 30], [20, 30, 40], [30, 40, 50]]
    for x in range(len(UTR_end_list_filtered)-2):
        if curr_strand == '+':
            first = UTR_end_list_filtered[x]
            middle = UTR_end_list_filtered[x+1]
            last = UTR_end_list_filtered[x+2]
        elif curr_strand == '-':
            first = UTR_end_list_filtered[x]
            middle = UTR_end_list_filtered[x+1]
            last = UTR_end_list_filtered[x+2]
        UTR_search_region.append([first,middle,last])
    return UTR_search_region

def UTR_end_list_filtering_backward(UTR_end_list):
    UTR_end_list = UTR_end_list[::-1]
    UTR_end_list_filtered = [UTR_end_list[x] for x in range(len(UTR_end_list)) if UTR_end_list[x] == UTR_end_list[0] or (UTR_end_list[x] == UTR_end_list[-1] and abs(UTR_end_list[-2] - UTR_end_list[-1]) + 1 >= 200) or abs(UTR_end_list[x+1] - UTR_end_list[x]) + 1 >= 200]
    UTR_end_list_filtered = UTR_end_list_filtered[::-1]

    #Reverse a list (if curr_stand is '-')
    if curr_strand == '+':
        UTR_end_list_filtered = UTR_end_list_filtered[::-1]
    elif curr_strand == '-':
        UTR_end_list_filtered = UTR_end_list_filtered[::-1]

    #Make triple dataset
    #UTR_end_list = [10,20,30,40,45]
    #UTR_start_site = 0
    #UTR_end_site = 50
    #UTR_search_region: [[0, 10, 20], [10, 20, 30], [20, 30, 40], [30, 40, 50]]
    for x in range(len(UTR_end_list_filtered)-2):
        if curr_strand == '+':
            first = UTR_end_list_filtered[x]
            middle = UTR_end_list_filtered[x+1]
            last = UTR_end_list_filtered[x+2]
        elif curr_strand == '-':
            first = UTR_end_list_filtered[x]
            middle = UTR_end_list_filtered[x+1]
            last = UTR_end_list_filtered[x+2]
        UTR_search_region.append([first,middle,last])
    return UTR_search_region


