#!/usr/bin/env python3

import os

def Merge_pA_site_infor(pA_site_bed_file, output_file):
    output = open(output_file, 'w')
    temp_dict = {}
    pA_dict = {}
    for line in open(pA_site_bed_file, 'r'):
        line = line.rstrip()
        fields = line.split("\t")
        bed_infor = fields[0:6]
        name = fields[3]
        pA_site = fields[8]
        temp_dict[name] = bed_infor
        if not name in pA_dict:
            pA_dict[name] = [pA_site]
        else:
            pA_dict[name].append(pA_site)

    for name in temp_dict.keys():
        bed_infor = temp_dict[name]
        pA_sites = '|'.join(pA_dict[name])
        print("\t".join(bed_infor),pA_sites, sep="\t", end="\n", file=output)
