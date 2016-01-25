#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def coverage_comparison_with_pA_site(All_Samples_curr_3UTR_coverages, UTR_start, UTR_end, curr_strand, weight_for_second_coverage, Coverage_pPAS_cutoff, pA_site, test_name):
    ###For each gene###
    #Parameter setting
    coverage_threshold = Coverage_pPAS_cutoff #Depth(Coverage) threshold in 5'end of last exon
    search_point_start = 200
    search_point_end = 200
    coverage_test_region = 100 #testing 0-100bp in 5'end of last exon
    least_3UTR_length = 500 #>=150bp 3'UTR length are needed

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
        #Prepare UTR search region
        UTR_start_site = UTR_start - 1
        UTR_end_site = UTR_end

        UTR_end_list = pA_site
        #TEST:
        #UTR_end_list = [10,20,30,40,45]
        #UTR_start_site = 0
        #UTR_end_site = 50
        #UTR_search_region: [[0, 10, 20], [10, 20, 30], [20, 30, 40], [30, 40, 50]]
        UTR_search_region = []
        if curr_strand == '+':
            if UTR_start_site < UTR_end_list[0]:
                UTR_end_list.insert(0,UTR_start_site)
            if UTR_end_list[-1] < UTR_end_site:
                if (UTR_end_site - UTR_end_list[-1]) <= 50: #Remove pA site near 3'end(Default: 50bp)
                    UTR_end_list[-1] = UTR_end_site
                else:
                    UTR_end_list.append(UTR_end_site)
        elif curr_strand == '-':
            if UTR_start_site < UTR_end_list[0]:
                if (UTR_end_list[0]- UTR_start_site) <= 50: #Remove pA site near 3'end(Default: 50bp)
                    UTR_end_list[0] = UTR_start_site
                else:
                    UTR_end_list.insert(0,UTR_start_site)
            if UTR_end_list[-1] < UTR_end_site:
                UTR_end_list.append(UTR_end_site)

        for x in range(len(UTR_end_list)-2):
            first = UTR_end_list[x]
            middle = UTR_end_list[x+1]
            last = UTR_end_list[x+2]
            UTR_search_region.append([first,middle,last])

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
