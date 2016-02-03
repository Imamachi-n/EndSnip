#!/usr/bin/env python3

import math
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

###STEP1###
def De_novo_UTR_search_region(UTR_start, UTR_end, curr_strand, search_region_distance):
    UTR_end_list = []
    UTR_search_region = []
    UTR_length = UTR_end - UTR_start + 1

    #Define the end sites of sub-UTR regions(chromosome sites)
    #UTR_start = 8023457
    #UTR_end = 8028691
    if curr_strand == '+':
        UTR_end_list = [x*search_region_distance + UTR_start for x in range(1, math.floor(UTR_length/search_region_distance))]
        UTR_end_list.insert(0,UTR_start)
        UTR_end_list.append(UTR_end)
    elif curr_strand == '-':
        UTR_end_list = [UTR_end - x*search_region_distance for x in range(1, math.floor(UTR_length/search_region_distance))]
        UTR_end_list.insert(0,UTR_end)
        UTR_end_list.append(UTR_start)

    #Define UTR search sub-region(Chromosome sites)
    for x in range(len(UTR_end_list)-2):
        if curr_strand == '+':
            first = UTR_end_list[x]
            middle = UTR_end_list[x+1]
            last = UTR_end_list[x+2]
        elif curr_strand == '-':
            first = UTR_end_list[x]
            middle = UTR_end_list[x+1]
            last = UTR_end_list[x+2]
        UTR_search_region.append([first,middle,last]) #Middle is meaningless.

    return UTR_search_region


###STEP2&3###
def Define_UTR_search_region(UTR_start_site, UTR_end_site, curr_strand, UTR_end_list):
    #Define 3'end site from pA_site and 3'end infromation from RefSeq(or the other transcriptome DB)
    UTR_search_region = []
    #if curr_strand == '+':
    #    if UTR_start_site < UTR_end_list[0]:
    #        UTR_end_list.insert(0,UTR_start_site)
    #    if UTR_end_list[-1] < UTR_end_site:
    #        if (UTR_end_site - UTR_end_list[-1]) <= 200: #Remove pA site near 3'end(Default: 50bp)
    #            UTR_end_list[-1] = UTR_end_site
    #        else:
    #            UTR_end_list.append(UTR_end_site)
    #elif curr_strand == '-':
    #    if UTR_start_site < UTR_end_list[0]:
    #        if (UTR_end_list[0] - UTR_start_site) <= 200: #Remove pA site near 3'end(Default: 50bp)
    #            UTR_end_list[0] = UTR_start_site
    #        else:
    #            UTR_end_list.insert(0,UTR_start_site)
    #    if UTR_end_list[-1] < UTR_end_site:
    #        UTR_end_list.append(UTR_end_site)

    #if curr_strand == '+':
    #    UTR_end_list = UTR_end_list[::-1]
    #elif curr_strand == '-':
    #    pass

    if curr_strand == '+':
        UTR_end_list.insert(0,UTR_start_site)
        UTR_end_list.append(UTR_end_site)
    elif curr_strand == '-':
        UTR_end_list.insert(0,UTR_end_site)
        UTR_end_list.append(UTR_start_site)

    print(UTR_end_list)
    #UTR_search_region = UTR_end_list_filtering_backward(UTR_end_list)
    UTR_search_region = Define_UTR_search_region_clustering(UTR_end_list)

    #Consider Read length(single-end 36bp) #TODO: 3'endのだらだらしたすそ部分の影響を緩和する。
    #UTR_search_region = [[8028691, 0, 8027503], [8028340, 0, 8023457]]
    #read_length = 36
    #if curr_strand == '+':
    #    UTR_search_region_re = [[(three_set[0] + read_length), three_set[1], three_set[2]] for three_set in UTR_search_region if not three_set[0] == UTR_search_region[0][0]]
    #    UTR_search_region_re.insert(0,UTR_search_region[0])
    #elif curr_strand == '-':
    #    UTR_search_region_re = [[(three_set[0] - read_length), three_set[1], three_set[2]] for three_set in UTR_search_region if not three_set[0] == UTR_search_region[0][0]]
    #    UTR_search_region_re.insert(0,UTR_search_region[0])

    return UTR_search_region

def Define_UTR_search_region_clustering(UTR_end_list):
    #UTR_end_list = [1, 201, 300, 600, 1000, 1500, 1600]
    #UTR_end_list = [89725044, 89725208, 89725774, 89726079, 89726202, 89726744, 89728093, 89728503, 89728544, 89731687]
    #UTR_end_list = [8028691, 8028417, 8028340, 8027503, 8027341, 8027157, 8023457]
    #[1 | 201, 300 | 600 | 1000 | 1500, 1600]
    #[ [1, 0, 600], [300, 0, 1000], [600, 0, 1500] ]
    test_list = []
    flg = 1
    for x in range(len(UTR_end_list)-1):  
        if x == 0:
            test_list.append(flg)
        elif abs(UTR_end_list[x+1] - UTR_end_list[x])+1 >= 200:
            test_list.append(flg)
        else:
            pass
        flg += 1
    test_list.insert(0,0)
    test_list.append(len(UTR_end_list)) #[0, 1, 3, 4, 5, 7]
    test_index = [[test_list[x],test_list[x+1]] for x in range(len(test_list)-1)] #[[0, 1], [1, 3], [3, 4], [4, 5], [5, 7]]
    #print(test_index)
    UTR_end_list_clustering = [UTR_end_list[test_index[x][0]:test_index[x][1]] for x in range(len(test_index))] #[[1], [201, 300], [600], [1000], [1500, 1600]]
    #print(UTR_end_list_clustering)
    UTR_search_region = [[UTR_end_list_clustering[x][-1],0,UTR_end_list_clustering[x+2][0]] for x in range(len(UTR_end_list_clustering)-2)] #[[1, 0, 600], [300, 0, 1000], [600, 0, 1500]]
    return UTR_search_region


###STEP: DIFF###
def Estimate_UTR_isoform_expression(curr_3UTR_curr_sample_bp_coverage, break_point_for_diff, curr_strand, chrom, test_wig_output_file1, test_wig_output_file2, flg_test):
    #Prepare index for each UTR region
    start_site = break_point_for_diff[0]
    UTR_segments = [[abs(break_point_for_diff[x]-start_site), abs(break_point_for_diff[x+1]-start_site)] for x in range(len(break_point_for_diff)-1)]
    #[[0, 274], [274, 1290], [1290, 5234]]
    if curr_strand == '+':
        segments_chrom_site = [[break_point_for_diff[x], break_point_for_diff[x+1]] for x in range(len(break_point_for_diff)-1)]
    elif curr_strand == '-':
        segments_chrom_site = [[break_point_for_diff[x+1], break_point_for_diff[x]] for x in range(len(break_point_for_diff)-1)]
    
    #multi-UTR coverages(median) #TODO: UTRの発現量(Coverage)の計算をMedianで行ってよいか要検討。
    multi_UTR_coverage = []
    for x in range(len(UTR_segments)):
        start_index = UTR_segments[x][0]
        end_index = UTR_segments[x][1]
        
        start_chrom_site = segments_chrom_site[x][0]
        end_chrom_site = segments_chrom_site[x][1]
        curr_UTR_bp_coverage = np.array(curr_3UTR_curr_sample_bp_coverage[start_index:end_index])
        curr_UTR_median_coverage = np.median(curr_UTR_bp_coverage)
        ###TEST: Wig file preparation
        if flg_test == 0: #siCTRL
            print(chrom, start_chrom_site, end_chrom_site, curr_UTR_median_coverage, sep="\t", end="\n", file=test_wig_output_file1)
        elif flg_test == 1: #siCFIm25
            print(chrom, start_chrom_site, end_chrom_site, curr_UTR_median_coverage, sep="\t", end="\n", file=test_wig_output_file2)
        multi_UTR_coverage.append(curr_UTR_median_coverage)

    #Estimate each UTR coverage
    Each_UTR_coverage = [multi_UTR_coverage[x] - multi_UTR_coverage[x+1] for x in range(len(multi_UTR_coverage)) if not x == (len(multi_UTR_coverage)-1)]
    Each_UTR_coverage.append(multi_UTR_coverage[-1]) #Required: [75.0, 39.0, 14.0, 2.0] => [36.0, 25.0, 12.0, 2.0]
    Each_UTR_coverage_sum = np.sum(np.array(Each_UTR_coverage))
    Each_UTR_coverage_percentage = np.array(Each_UTR_coverage) / Each_UTR_coverage_sum #[0.47999999999999998, 0.33333333333333331, 0.16, 0.026666666666666668]

    return [multi_UTR_coverage, Each_UTR_coverage, list(Each_UTR_coverage_percentage)]


def Estimate_PDUI_score(Each_UTR_coverage, Each_UTR_coverage_percentage, num_group_1, num_group_2):
    #Prepare Mean Each_UTR_coverage_percentage for two conditions
    Each_UTR_coverage_percentage_sample1 = np.array(Each_UTR_coverage_percentage[:num_group_1])
    Each_UTR_coverage_percentage_sample1_mean = np.mean(Each_UTR_coverage_percentage_sample1, axis=0)
    Each_UTR_coverage_percentage_sample2 = np.array(Each_UTR_coverage_percentage[num_group_1:])
    Each_UTR_coverage_percentage_sample2_mean = np.mean(Each_UTR_coverage_percentage_sample2, axis=0)

    #The difference of Mean Each_UTR_coverage_percentage between two conditions
    Each_UTR_coverage_percentage_sub = np.array(Each_UTR_coverage_percentage_sample2_mean) - np.array(Each_UTR_coverage_percentage_sample1_mean)

    ##Prepare Mean Each_UTR_coverage for two conditions
    #Each_UTR_coverage_sample1 = np.array(Each_UTR_coverage[:num_group_1])
    #Each_UTR_coverage_sample1_mean = np.mean(Each_UTR_coverage_sample1, axis=0)
    #Each_UTR_coverage_sample2 = np.array(Each_UTR_coverage[num_group_1:])
    #Each_UTR_coverage_sample2_mean = np.mean(Each_UTR_coverage_sample2, axis=0)

    #Define short-UTR/long-UTR isoform sets
    first_UTR_status = Each_UTR_coverage_percentage_sub[0]
    diff_index = 0
    short_UTR_coverage_sample1_mean = 0
    short_UTR_coverage_sample2_mean = 0
    for x in range(len(Each_UTR_coverage_percentage_sub)):
        if first_UTR_status < 0:
            if Each_UTR_coverage_percentage_sub[x] < 0:
                #short_UTR_coverage_sample1_mean += Each_UTR_coverage_sample1_mean[x]
                #short_UTR_coverage_sample2_mean += Each_UTR_coverage_sample2_mean[x]
                diff_index += 1
            else:
                break
        elif first_UTR_status >= 0:
            if Each_UTR_coverage_percentage_sub[x] >= 0:
                #short_UTR_coverage_sample1_mean += Each_UTR_coverage_sample1_mean[x]
                #short_UTR_coverage_sample2_mean += Each_UTR_coverage_sample2_mean[x]
                diff_index += 1
            else:
                break

    #Estimate PDUI score for each sample
    line_write = []
    PDUI_all = []
    short_UTR_coverage_all = []
    for x in range((num_group_1 + num_group_2)):
        curr_sample_each_UTR_coverage = Each_UTR_coverage[x]
        curr_sample_each_UTR_coverage_short = np.sum(np.array(curr_sample_each_UTR_coverage[:diff_index]))
        curr_sample_each_UTR_coverage_long = np.sum(np.array(curr_sample_each_UTR_coverage[diff_index:]))
        curr_sample_each_UTR_coverage_all = np.sum(np.array(curr_sample_each_UTR_coverage))
        PDUI_sample = curr_sample_each_UTR_coverage_long / curr_sample_each_UTR_coverage_all
        short_UTR_coverage_all.append(curr_sample_each_UTR_coverage_short)
        PDUI_all.append(PDUI_sample)
        line_write.extend([curr_sample_each_UTR_coverage_long, curr_sample_each_UTR_coverage_short, PDUI_sample])

    #Estimate mean PDUI and diff-PDUI
    short_UTR_coverage_sample1_mean = np.mean(np.array(short_UTR_coverage_all[:num_group_1]))
    short_UTR_coverage_sample2_mean = np.mean(np.array(short_UTR_coverage_all[num_group_1:]))
    PDUI_sample1_mean = np.mean(np.array(PDUI_all[:num_group_1]))
    PDUI_sample2_mean = np.mean(np.array(PDUI_all[num_group_1:]))
    PDUI_diff = PDUI_sample2_mean - PDUI_sample1_mean
    Fold_change = PDUI_sample1_mean / PDUI_sample2_mean
    Fold_chamge2 = short_UTR_coverage_sample2_mean / short_UTR_coverage_sample1_mean
    line_write.extend([PDUI_sample1_mean, PDUI_sample2_mean, PDUI_diff, Fold_change, Fold_chamge2])

    #Casting
    line_write = map(str, line_write)

    return line_write

    #short_UTR_sum_coverage_sample1 = 0
    #short_UTR_sum_coverage_sample2 = 0

    ##Estimate all UTR-isoforms coverage
    #all_UTR_sum_coverage_sample1 = np.sum(np.array(Each_UTR_coverage[:num_group_1]), axis=0)
    #all_UTR_sum_coverage_sample2 = np.sum(np.array(Each_UTR_coverage[num_group_1:]), axis=0)
    ##all_UTR_sum_coverage_sample1 = np.sum(np.array(Each_UTR_coverage[0]))
    ##all_UTR_sum_coverage_sample2 = np.sum(np.array(Each_UTR_coverage[1]))

    #Each_UTR_coverage_sub = np.array(all_UTR_sum_coverage_sample2) - np.array(all_UTR_sum_coverage_sample1)
    #first_UTR = Each_UTR_coverage_percentage_sub[0]

    #diff_index = 0
    #for x in range(len(Each_UTR_coverage_sub)):
    #    if first_UTR < 0:
    #        if Each_UTR_coverage_sub[x] < 0:
    #            short_UTR_sum_coverage_sample1 += all_UTR_sum_coverage_sample1[0][x]
    #            short_UTR_sum_coverage_sample2 += all_UTR_sum_coverage_sample1[1][x]
    #        else:
    #            break
    #    elif first_UTR >= 0:
    #        if Each_UTR_coverage_sub[x] >= 0:
    #            short_UTR_sum_coverage_sample1 += Each_UTR_coverage[0][x]
    #            short_UTR_sum_coverage_sample2 += Each_UTR_coverage[1][x]
    #        else:
    #            break

    #Estimate PDUI score
    #PDUI_sample1 = (all_UTR_sum_coverage_sample1 - short_UTR_sum_coverage_sample1) / all_UTR_sum_coverage_sample1
    #PDUI_sample2 = (all_UTR_sum_coverage_sample2 - short_UTR_sum_coverage_sample2) / all_UTR_sum_coverage_sample2

    #PDUI_diff = PDUI_sample2 - PDUI_sample1
    #Fold_change = PDUI_sample1 / PDUI_sample2
    #Fold_change2 = short_UTR_sum_coverage_sample2 / short_UTR_sum_coverage_sample1

    #print(PDUI_sample1, PDUI_sample2, PDUI_diff, Fold_change, Fold_change2)



###########################################
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
            Evaluation_result = Check_break_point(curr_UTR_search_coverage, min_MSE_index, search_point_start, search_point_end, curr_strand)

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

def Check_break_point(curr_UTR_search_coverage, min_MSE_index, search_point_start, search_point_end, curr_strand):
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

    if div_median >= 1.5: #TODO: LongUTRとShortUTRの差異に関して、Criteriaをどうするか決める。
        return True
    else:
        return False

###########################################
