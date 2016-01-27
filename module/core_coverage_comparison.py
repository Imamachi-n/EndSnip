#!/usr/bin/env python3

import math
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt


def Define_UTR_search_region(UTR_start_site, UTR_end_site, curr_strand, UTR_end_list):
    #Define 3'end site from pA_site and 3'end infromation from RefSeq(or the other transcriptome DB)
    UTR_search_region = []
    if curr_strand == '+':
        if UTR_start_site < UTR_end_list[0]:
            UTR_end_list.insert(0,UTR_start_site)
        if UTR_end_list[-1] < UTR_end_site:
            if (UTR_end_site - UTR_end_list[-1]) <= 200: #Remove pA site near 3'end(Default: 50bp)
                UTR_end_list[-1] = UTR_end_site
            else:
                UTR_end_list.append(UTR_end_site)
    elif curr_strand == '-':
        if UTR_start_site < UTR_end_list[0]:
            if (UTR_end_list[0] - UTR_start_site) <= 200: #Remove pA site near 3'end(Default: 50bp)
                UTR_end_list[0] = UTR_start_site
            else:
                UTR_end_list.insert(0,UTR_start_site)
        if UTR_end_list[-1] < UTR_end_site:
            UTR_end_list.append(UTR_end_site)

    if curr_strand == '+':
        UTR_end_list = UTR_end_list[::-1]
    elif curr_strand == '-':
        pass

    #Remove pA_site near the other pA_site(Default: >=150bp)
    ###TEST(Strand:-): UTR_end_list = [1,200,400,450,600,800,850,890,1000]
    ###*      *      -   *        *    -   -   -    *  check
    ###|      |      |   |        |    |   |   |    |
    ###1     200    400 450      600  800 850 890 1000 chrom_site

    ###TEST(Strand:+): UTR_end_list = [1000,800,610,550,500,300,1]
    ###  *       *      -   -   *      *        * check
    ###  |       |      |   |   |      |        |
    ###1000     800    610 550 500    300       1 chrom_site
    UTR_end_list_filtered = [UTR_end_list[x] for x in range(len(UTR_end_list)) if UTR_end_list[x] == UTR_end_list[-1] or abs(UTR_end_list[x+1] - UTR_end_list[x] + 1) >= 150]

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


def De_novo_UTR_search_region(UTR_start, UTR_end, curr_strand, search_region_distance):
    UTR_end_list = []
    UTR_search_region = []
    UTR_length = UTR_end - UTR_start + 1

    for x in range(1, math.floor(UTR_length)+1):
        UTR_end_list = [x*search_region_distance + UTR_start for x in range(1, math.floor(UTR_length/search_region_distance))]
        UTR_end_list.insert(0,UTR_start)
        UTR_end_list.append(UTR_end)

    if curr_strand == '-':
        UTR_end_list = UTR_end_list[::-1]

    for x in range(len(UTR_end_list)-2):
        if curr_strand == '+':
            first = UTR_end_list[x]
            middle = UTR_end_list[x+1]
            last = UTR_end_list[x+2]
        elif curr_strand == '-':
            first = UTR_end_list[x]
            middle = UTR_end_list[x+1]
            last = UTR_end_list[x+2]
        UTR_search_region.append([first,middle,last])

    return UTR_search_region


def Estimate_break_point(curr_3UTR_curr_sample_bp_coverage, UTR_search_region, UTR_start_site, UTR_end_site, curr_strand, search_point_start, search_point_end, least_search_region_coverage, test_name, flg):
    #Initiate break point candidate list
    break_point_candidates = [] #Chrom_site

    for curr_UTR_search_region in UTR_search_region:
        #Initialize variance and estimated 3'UTR abundance lists
        #Initialize break point infor list
        curr_mean_variance_list = []
        curr_estimated_3UTR_abundance_list = []

        #Define UTR coverage for each UTR search region
        curr_UTR_search_coverage = []
        search_region_length = 0
        if curr_strand == '+':
            start_site = curr_UTR_search_region[0] - UTR_start_site
            end_site = curr_UTR_search_region[2] - UTR_start_site + 1 + 1 #length(end-start)+1, slicing
            curr_UTR_search_coverage = curr_3UTR_curr_sample_bp_coverage[start_site:end_site]
            search_region_length = end_site - start_site
            print(start_site,end_site,search_region_length)
        elif curr_strand == '-':
            start_site = abs(curr_UTR_search_region[0] - UTR_end_site)
            end_site = abs(curr_UTR_search_region[2] - UTR_end_site) + 1 + 1 #length(end-start)+1, slicing
            curr_UTR_search_coverage = curr_3UTR_curr_sample_bp_coverage[start_site:end_site]
            search_region_length = end_site - start_site
            print(start_site,end_site,search_region_length)

        #TODO: TEST: search region coverage variance
        #curr_UTR_search_coverage_mean = np.mean(np.array(curr_UTR_search_coverage))
        #curr_UTR_search_coverage_variance = (curr_UTR_search_coverage - curr_UTR_search_coverage_mean)**2
        #curr_UTR_search_coverage_variance_mean = np.mean(curr_UTR_search_coverage_variance)

        #Define coverage variance search region #TEST: 1-700bp(Index: 0-699)/700bp length => Index: 100:600(101-600bp)
        start_site_for_coverage_variance_search = search_point_start #Index: 100
        end_site_for_coverage_variance_search = search_region_length - search_point_end #Index: 600

        #Initiate a list of results for each sample
        curr_region_result = [ [], [], [] ] # Mean_squared_error | Long_UTR_coverage | Short_UTR_coverage

        #Check mean squared error(variance) in search region
        for curr_point in range(start_site_for_coverage_variance_search, end_site_for_coverage_variance_search):
            curr_break_point = curr_point #Current base-pair on 3UTR region for checking

            #Testing coverage variance for each base-pair(bp)
            Mean_squared_error = 0
            Long_UTR_coverage = 0
            Short_UTR_coverage = 0
            Mean_squared_error, Long_UTR_coverage, Short_UTR_coverage = Estimate_variance(curr_UTR_search_coverage, curr_break_point)
            #if Mean_squared_error is None:
            #    continue
            #else:
            curr_region_result[0].append(Mean_squared_error)
            curr_region_result[1].append(Long_UTR_coverage)
            curr_region_result[2].append(Short_UTR_coverage)

        #TODO: TEST: Mean_squared_error in 3'UTR region for PTEN, ELAVL1
        name = 'CTRL'
        if flg == 1:
            name = 'CFIm25KD'
        plt.plot(curr_region_result[0])
        #plt.show()
        filename = "data/output_variance_" + test_name + '_' + name + '_' + str(curr_UTR_search_region[0]) + "-" + str(curr_UTR_search_region[2]) + ".png"
        plt.savefig(filename)
        plt.close()

        #plt.plot(curr_UTR_search_coverage_variance)
        #filename = "data/output_variance_mean_" + test_name + '_' + name + '_' + str(curr_UTR_search_region[0]) + "-" + str(curr_UTR_search_region[2]) + ".png"
        #plt.savefig(filename)
        #plt.close()

        #Identify index of min mean squared error in curr_region_result[0]
        min_MSE_index = curr_region_result[0].index(min(curr_region_result[0])) #Break point
        break_point_chrom_site = 0
        if curr_strand == '+':
            break_point_chrom_site = curr_UTR_search_region[0] + min_MSE_index + search_point_start
        elif curr_strand == '-':
            break_point_chrom_site = curr_UTR_search_region[0] - min_MSE_index - search_point_start
        break_point_long_UTR_coverage = curr_region_result[1][min_MSE_index]
        break_point_short_UTR_coverage = curr_region_result[2][min_MSE_index]

        #CHECK: least search region coverage(Default: >=5 coverage)
        #CHECK: Short UTR coverage > Long UTR coverage
        #CHECK: 
        if break_point_short_UTR_coverage >= least_search_region_coverage and break_point_short_UTR_coverage > break_point_long_UTR_coverage:
            #Evaluate difference between short and long 3UTR(Default: At least 2-fold changed)
            Evaluation_result = Evaluate_break_point(curr_UTR_search_coverage, min_MSE_index, search_point_start, search_point_end, curr_strand)

            if Evaluation_result:
                res = break_point_short_UTR_coverage - break_point_long_UTR_coverage
                dev = break_point_short_UTR_coverage / break_point_long_UTR_coverage
                print(break_point_chrom_site, break_point_short_UTR_coverage, break_point_long_UTR_coverage, res, dev)

                #Append Chrom sites of break point candidates
                break_point_candidates.append(break_point_chrom_site)
            else:
                print("NG!!")
        else:
            print("NG!!")

    return break_point_candidates
            
def Estimate_variance(curr_UTR_search_coverage, curr_break_point):
    search_coverage_long = np.array(curr_UTR_search_coverage[curr_break_point:])
    search_coverage_short = np.array(curr_UTR_search_coverage[0:curr_break_point])

    #Calculate mean of coverage(Long/Short 3UTR)
    Long_UTR_coverage = np.mean(search_coverage_long)
    Short_UTR_coverage = np.mean(search_coverage_short)

    #if Short_UTR_coverage > Long_UTR_coverage:
    search_coverage_long_residual = search_coverage_long - Long_UTR_coverage
    search_coverage_short_residual = search_coverage_short - Short_UTR_coverage

    coverage_residual = search_coverage_long_residual
    coverage_residual = np.append(coverage_residual, search_coverage_short_residual)
    Mean_squared_error = np.mean(coverage_residual**2)
    
    return Mean_squared_error, Long_UTR_coverage, Short_UTR_coverage
    #else:
    #    return None, None, None #Null objects

def Evaluate_break_point(curr_UTR_search_coverage, min_MSE_index, search_point_start, search_point_end, curr_strand):
    best_break_point = min_MSE_index + search_point_start + 1 #Slicing
    #print(best_break_point)

    #Extract search coverages in short/long 3UTR region
    search_coverage_long = np.array(curr_UTR_search_coverage[best_break_point:])
    search_coverage_short = np.array(curr_UTR_search_coverage[0:best_break_point])

    #print(search_coverage_short)
    #Wilcoxon test??
    #pvalue_test1 = scipy.stats.ranksums(search_coverage_long, search_coverage_short)
    #pvalue_test2 = scipy.stats.mannwhitneyu(search_coverage_long, search_coverage_short)
    #print(pvalue_test2)

    #Calculate mean of coverage(Long/Short 3UTR)
    #search_coverage_mean_long = np.mean(search_coverage_long)
    #search_coverage_mean_short = np.mean(search_coverage_short)
    search_coverage_median_long = np.median(search_coverage_long)
    search_coverage_median_short = np.median(search_coverage_short)

    #div_mean = search_coverage_mean_short / search_coverage_mean_long
    div_median = search_coverage_median_short / search_coverage_median_long
    print(div_median)

    #Residual
    #search_coverage_residual_long = (search_coverage_long - search_coverage_mean_long)
    #search_coverage_residual_short = (search_coverage_short - search_coverage_mean_short)

    if div_median >= 2: #TODO: LongUTRとShortUTRの差異に関して、Criteriaをどうするか決める。
        return True
    else:
        return False


###MAIN###
def de_novo_coverage_comparison_with_windows(curr_3UTR_all_samples_bp_coverage, curr_3UTR_all_samples_bp_chrom_site, UTR_start, UTR_end, curr_strand, weight_for_second_coverage, Coverage_pPAS_cutoff, pA_site, test_name):
    ###For each gene###
    #Parameter setting
    coverage_threshold = Coverage_pPAS_cutoff #Depth(Coverage) threshold in 5'end of last exon
    search_point_start = 100
    search_point_end = 100
    coverage_test_region = 100 #testing 0-100bp in 5'end of last exon
    least_3UTR_length = 500 #>=150bp 3'UTR length are needed
    least_search_region_coverage = Coverage_pPAS_cutoff #>=5 coverage are needed in search region
    search_region_distance = 200

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

        #Estimate_each_sample
        flg = 0 #TODO: TEST:
        for curr_3UTR_curr_sample_bp_coverage in curr_3UTR_all_samples_bp_coverage:
            break_point_candidates = Estimate_break_point(curr_3UTR_curr_sample_bp_coverage, UTR_search_region, UTR_start_site, UTR_end_site, curr_strand, search_point_start, search_point_end, least_search_region_coverage, test_name, flg)

            if curr_strand == '-':
                break_point_candidates = break_point_candidates[::-1]
            print(break_point_candidates)

            #ReEstimate break points
            updated_UTR_search_region = Define_UTR_search_region(UTR_start_site, UTR_end_site, curr_strand, break_point_candidates)
            print(updated_UTR_search_region)
            updated_break_point_candidates = Estimate_break_point(curr_3UTR_curr_sample_bp_coverage, updated_UTR_search_region, UTR_start_site, UTR_end_site, curr_strand, search_point_start, search_point_end, least_search_region_coverage, test_name, flg)

            flg = 1 #TODO: TEST:


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
