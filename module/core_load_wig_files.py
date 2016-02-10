#!/usr/bin/env python3

from module.utils_others import now_time
from bisect import bisect
import numpy as np

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
        #pA_site = fields[6]
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
        
        region_start = int(region_start) + 1 #0-base => 1-base
        region_end = int(region_end)

        if (region_end - region_start) >= 500: #Min 3UTR length(Default: 500bp)
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
            #Load wig file
            chrom_name = fields[0]
            region_start = int(fields[1])
            region_end = int(fields[2])
            read_depth = int(float(fields[-1]))
            
            #Initialize coverage data in each chromosome
            if chrom_name not in curr_sample_All_chroms_coverage_dict:
                curr_sample_All_chroms_coverage_dict[chrom_name] = [[0],[0]] #[[region_site], [depth]]

            #Add coverage data in each region on each chromosome
            if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]: #if gap region exists
                curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start) #Region end => Region start #1-based
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
            region_start = curr_3UTR_structure[1] #3'UTR region start
            region_end = curr_3UTR_structure[2] #3'UTR region end

            #Call current chromosome from dictionary
            if curr_chr in curr_sample_All_chroms_coverage_dict.keys():
                #Region and Depth for current chromosome
                curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr]

                #TEST:
                #Raw_data
                #chrom_site   = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
                #chr_coverage = [1, 4, 4, 0, 0, 0, 1, 2, 2, 0]
                #Input_data
                #chrom_site   = [0, 1, 3, 6, 7, 9, 10]
                #chr_coverage = [0, 1, 4, 0, 1, 2, 0]
                #from bisect import bisect
                #curr_chr_coverage = [[0, 1, 3, 6, 7, 9, 10], [0, 1, 4, 0, 1, 2, 0]]

                #from bisect import bisect
                #curr_chr_coverage = [[1,10,20,30,40,50,60], [1,10,10,10,30,30,50]] #[[Chrom_site],[Coverage]]
                #NO1:
                #region_start = 1
                #region_end = 60
                #NO2:
                #region_start = 1
                #region_end = 55
                #NO3:
                #region_start = 5
                #region_end = 60
                #NO4:
                #region_start = 5
                #region_end = 55
                #NO5:
                #region_start = 70
                #region_end = 90
                left_region_index = bisect(curr_chr_coverage[0], region_start) #Insertion site(index) of region start
                right_region_index = bisect(curr_chr_coverage[0], region_end) #Insertion site(index) of region end

                extracted_3UTR_region = []
                extracted_coverage = []

                #In the case of 0 coverage,
                if left_region_index == right_region_index:
                    extracted_3UTR_region = [region_start, region_end]
                    extracted_coverage = [0, 0]
                elif int(curr_chr_coverage[0][left_region_index-1]) == int(region_start) and int(curr_chr_coverage[0][right_region_index-1]) == int(region_end):
                    #print("1")
                    #List of 3UTR region
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index-1:right_region_index]
                    #List of depth(coverage) in 3'UTR region
                    extracted_coverage = curr_chr_coverage[1][left_region_index-1:right_region_index]
                elif int(curr_chr_coverage[0][left_region_index-1]) == int(region_start):
                    #print("2")
                    #List of 3UTR region
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index-1:right_region_index]
                    extracted_3UTR_region.append(region_end)
                    #List of depth(coverage) in 3'UTR region
                    extracted_coverage = curr_chr_coverage[1][left_region_index-1:right_region_index+1]
                elif int(curr_chr_coverage[0][right_region_index-1]) == int(region_end):
                    #print("3")
                    #List of 3UTR region
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0,region_start)
                    #List of depth(coverage) in 3'UTR region
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index]
                    extracted_coverage.insert(0,curr_chr_coverage[1][left_region_index])
                else:
                    #print("4")
                    #List of 3UTR region
                    extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                    extracted_3UTR_region.insert(0,region_start)
                    extracted_3UTR_region.append(region_end)
                    #List of depth(coverage) in 3'UTR region
                    extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                    extracted_coverage.insert(0,curr_chr_coverage[1][left_region_index])

                '''
                #Example
                #Index:        0  1  2  3  4  5  6
                #chrom_site = [1,10,20,30,40,50,60]
                #coverage   = [1,10,10,10,30,30,50]
                #1bp => 1
                #2-10bp => 10
                #11-20bp => 10
                #21-30bp => 10
                #31-40bp => 30
                #41-50bp => 30
                #51-60bp => 50

                ###1-60[1,7]###
                #bisect(chrom_site,1) => 1 => 0/0
                #bisect(chrom_site,60) => 7 => 7/7
                #chrom_site: [1,10,20,30,40,50,60] => [1,10,20,30,40,50,60]
                #coverage:   [1,10,10,10,30,30,50] => [1,10,10,10,30,30,50]

                ###1-55[1,6]###
                #bisect(chrom_site,1) => 1 => 0/0
                #bisect(chrom_site,55) => 6 => 6/7
                #chrom_site: [1,10,20,30,40,50]    => [1,10,20,30,40,50,"55"]
                #coverage:   [1,10,10,10,30,30,50] => [1,10,10,10,30,30, 50]

                ###5-60[1,7]###
                #bisect(chrom_site,5) => 1 => 1/1
                #bisect(chrom_site,60) => 7 => 7/7
                #chrom_site: [10,20,30,40,50,60] => [ "5",10,20,30,40,50,60]
                #coverage:   [10,10,10,30,30,50] => ["10",10,10,10,30,30,50]

                ###5-55[1,7]###
                #bisect(chrom_site,5) => 1 => 1/1
                #bisect(chrom_site,55) => 6 => 6/7
                #chrom_site: [10,20,30,40,50]    => [ "5",10,20,30,40,50,"55"]
                #coverage:   [10,10,10,30,30,50] => ["10",10,10,10,30,30, 50 ]
                '''

                #List of depth(coverage) in 3'UTR region
                #extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                
                #List of 3UTR region
                #extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                #extracted_3UTR_region.insert(0, region_start)
                #extracted_3UTR_region.append(region_end)

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
    UTR_length = extracted_3UTR_region[-1] - extracted_3UTR_region[0] + 1
    bp_coverage = np.zeros(UTR_length)
    bp_chrom_site = list(range(extracted_3UTR_region[0], extracted_3UTR_region[-1] + 1))
    print(len(bp_coverage))
    print(len(bp_chrom_site))
    print(bp_chrom_site[0],bp_chrom_site[-1])

    #First bp
    bp_coverage[0] = extracted_coverage[0]

    #Start site in 3UTR region on chromosome
    relative_start = extracted_3UTR_region[0]

    #Make a list of 1bp-resolution coverage in 3UTR region from 5'end of last exon (0-base)
    for i in range(len(extracted_coverage)-1):
        curr_region_start = extracted_3UTR_region[i] - relative_start + 1 #start site in its transcript
        curr_region_end = extracted_3UTR_region[i+1] - relative_start + 1 #end site in its transcript
        bp_coverage[curr_region_start:curr_region_end] = extracted_coverage[i] #List of coverage
        
        '''
        #Example
        #chrom_site: [5,10,20,30,40,55]
        #coverage:   [1,10,10,10,30,30]

        #First bp
        #chrom_site: [5] (1bp)
        #coverage:   [1]

        #Next bps
        #chrom_site: [5,10] =(-5)> [0,5] =(+1)> [1:6] (2-6bp)
        #coverage: 10
        ...
        #chrom_site: [10,20] =(-5)> [5,15] =(+1)> [6:16] (7-16bp)
        #coverage: 10
        '''

    #If stand is minus '-', reverse 1bp-resolution coverage list
    if strand_info == '-':
        bp_coverage = bp_coverage[::-1]

    #print((bp_coverage[0:10]))
    #print((bp_chrom_site))
    bp_coverage = list(bp_coverage)
    return [bp_coverage,bp_chrom_site]
