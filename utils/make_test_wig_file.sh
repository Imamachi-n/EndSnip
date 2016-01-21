#!/bin/bash
#python3 make_test_wig_file.py ../refGene.txt PTEN,ELAVL1 > test_genes.bed
bedtools intersect -a accepted_hits.bam.wig -b ../test_genes.bed -wa > accepted_hits_PTEN_ELAVL1.bam.wig