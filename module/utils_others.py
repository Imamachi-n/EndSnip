#!/usr/bin/env python3

def Input_config_file(config_file):
    config_dict = {}
    for line in open(config_file, 'r'):
        line = line.rstrip()
        if line == '' or line == '#':
            continue
        fields = line.split('=')
        key_name = fields[0]
        infor = fields[1]
        config_dict[key_name] = infor
    return config_dict

from datetime import datetime
#Print out comment with now time
def now_time(comment):
    now = datetime.now()
    nowtime = "{0:%Y-%m-%d %H:%M:%S}".format(now)
    print ('[' + nowtime + ']' + ' ' + comment)



#Parse configure file
def parse_cfgfile(cfg_file):
    Group1_Tophat_aligned_file = ''
    Group2_Tophat_aligned_file = ''
    output_directory = ''
    Annotated_3UTR_file = ''
    PolyA_site_infor = ''
    Output_result_file = ''

    Num_least_in_group1_local = ''
    Num_least_in_group2_local = ''
    Coverage_cutoff_local = ''
    FDR_cutoff_local = ''
    Fold_change_cutoff_local = ''
    PDUI_cutoff_local = ''
    Coverage_pPAS_cutoff_local = ''

    for line in open(cfg_file, 'r'):
        if line == '\n' or line == '#':
            continue
        line = line.rstrip()
        line = line.replace(' ','')
        command = line.split("=")

        #File/Directory
        if command[0] == 'Group1_Tophat_aligned_Wig':
            Group1_Tophat_aligned_file = command[1].split(',')
        if command[0] == 'Group2_Tophat_aligned_Wig':
            Group2_Tophat_aligned_file = command[1].split(',')
        if command[0] == 'Output_directory':
            output_directory = command[1]
            if output_directory[-1] != '/':
                output_directory += '/'
        if command[0] == 'Annotated_3UTR':
            Annotated_3UTR_file = command[1]
        if command[0] == 'PolyA_site_infor':
            PolyA_site_infor = command[1]
        if command[0] == 'Output_result_file':
            Output_result_file = command[1]

        #Parameter
        if command[0] == 'Num_least_in_group1':
            Num_least_in_group1_local = command[1]
        if command[0] == 'Num_least_in_group2':
            Num_least_in_group2_local = command[1]
        if command[0] == 'Coverage_cutoff':
            Coverage_cutoff_local = command[1]
        if command[0] == 'FDR_cutoff':
            FDR_cutoff_local = command[1]
        if command[0] == 'PDUI_cutoff':
            Fold_change_cutoff_local = command[1]
        if command[0] == 'Fold_change_cutoff':
            PDUI_cutoff_local = command[1]
        if command[0] == 'Coverage_pPAS_cutoff':
            Coverage_pPAS_cutoff_local = command[1]

    #Error handing
    if Group1_Tophat_aligned_file == '':
        sys.exit("ERROR: No Tophat aligned BAM file for group 1...")
    if Group2_Tophat_aligned_file == '':
        sys.exit("ERROR: No Tophat aligned BAM file for group 2...")
    if output_directory == '':
        sys.exit("ERROR: No output directory...")
    if Annotated_3UTR_file == '':
        sys.exit("ERROR: No annotated 3'UTR file...")
    if PolyA_site_infor == '':
        sys.exit("ERROR: No polyA site information file...")
    if Output_result_file == '':
        sys.exit("ERROR: No result file name...")

    return Group1_Tophat_aligned_file, Group2_Tophat_aligned_file, output_directory, Annotated_3UTR_file, PolyA_site_infor, Output_result_file, Num_least_in_group1_local, Num_least_in_group2_local, Coverage_cutoff_local, FDR_cutoff_local, Fold_change_cutoff_local, PDUI_cutoff_local, Coverage_pPAS_cutoff_local
