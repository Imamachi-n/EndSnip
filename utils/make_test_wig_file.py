#!/usr/bin/env python3
'''
Usage:
  make_test_wig_file.py [refGene_file(UCSC dataset)] [comma-separated Gene list(e.g. PTEN,ELAVL1)] > Selected_wig_file(Output_file)
'''

import sys

refGene_file = open(sys.argv[1], 'r') #../refGene.txt
gene_list = sys.argv[2].split(',') #PTEN,ELAVL1

ref_dict = {}

for line in refGene_file:
    line = line.rstrip()
    fields = line.split("\t")
    name = fields[12]
    if name in gene_list and not name in ref_dict:
        chrom = fields[2]
        st = int(fields[4]) - 1000
        ed = int(fields[5]) + 1000
        print(chrom,str(st),str(ed),sep="\t",end="\n")
        ref_dict[name] = 1
