#!/usr/bin/env python3

import sys
import re

input_file = open(sys.argv[1], 'r')
output_file = open(sys.argv[2], 'w')

if input_file == '' or output_file == '':
    sys.exit("ERROR: Indicate input/output files...")

for line in input_file:
    data = line.rstrip().split("\t")
    feature = data[2]
    infor = data[8].split("; ")
    if re.match('^transcript_id',infor[2]):
        name = infor[2].replace('transcript_id "','')
        name = name.replace('"','')
    elif re.match('^transcript_id',infor[3]):
        name = infor[3].replace('transcript_id "','')
        name = name.replace('"','')
    else:
        print("ERROR: "+infor[0])
        name = ''
    if re.match('^NM_',name) and feature == 'exon':
        continue
    elif re.match('^NM_',name) and feature == 'CDS':
        data[2] = 'exon'
    print("\t".join(data), end="\n",file=output_file)
