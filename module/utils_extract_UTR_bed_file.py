#!/usr/bin/env python3

import re

def Extract_gene_symbol_map_kfXref_file(gene_gtf_file):
    ref_dict = {}
    for line in open(gene_gtf_file, 'r'):
        line = line.rstrip()
        fields = line.split("\t")
        infor = fields[8].split("; ")
        gene_symbol = ''
        trx_id = ''
        for x in infor:
            if re.match('^gene_name', x):
                gene_symbol = x.replace('gene_name "', '')
                gene_symbol = gene_symbol.replace('"', '')
            if re.match('^transcript_id', x):
                trx_id = x.replace('transcript_id "', '')
                trx_id = trx_id.replace('"', '')
        if gene_symbol != '' and trx_id != '':
            ref_dict[trx_id] = gene_symbol
    return ref_dict        

##TODO: multi-3'UTRのケースに対応させる。
def Extract_3UTR_from_bed(gene_bed_file, gene_symbol_map_kfXref_file, output_utr_file):
    output_write = open(output_utr_file, 'w') #Output file
    
    raw_utr_dict = {} #UTR information dictionary
    refseq_trx_gene_symbol_dict = {} #RefSeq id => Gene Symbol dictionary
    num_line = 0 #Counter

    #Prepare RefSeq id => Gene Symbol dictionary
    for line in open(gene_symbol_map_kfXref_file, 'r'):
        if line == "\n" or line == '#':
            continue
        fields = line.strip().split("\t")
        gene_symbol = fields[1]
        refseq_trx_id = fields[0]
        refseq_trx_gene_symbol_dict[refseq_trx_id] = gene_symbol

    #for each transcript
    scanned_3UTR_list = {}
    num_saved = 0
    for line in open(gene_bed_file, 'r'):
        fields = line.rstrip().split("\t")
        refseq_id = fields[3] #refseq_id
        if '_' not in fields[0]: #chr1, chr2...
            if not refseq_id in refseq_trx_gene_symbol_dict:
                gene_symbol = "NA"
            else:
                gene_symbol = refseq_trx_gene_symbol_dict[refseq_id]

            UTR_id = [refseq_id, gene_symbol, fields[0], fields[5]] #refseq_id - gene_symbol - chrom - strand
            UTR_id_new = '|'.join(UTR_id)

            curr_strand = fields[5]
            gene_start = int(fields[1]) #Transcript start site
            UTR_start = ''
            UTR_end = ''
            this_UTR = ''
            UTR_end_new = ''
            if curr_strand == '+': #Strand: +
                UTR_start = str(gene_start + int(fields[-1].strip(',').split(',')[-1])) #UTR start(0-base): Trx start site + UTR start position in trx
                UTR_end = str(int(UTR_start) + int(fields[-2].strip(',').split(',')[-1])) #UTR end: UTR start + last exon length
                this_UTR = fields[0] + UTR_start + curr_strand #Chrom - UTR start site - stand(+/-)
                UTR_end_new = UTR_end
            elif curr_strand == '-': #Strand: -
                UTR_start = str(gene_start + int(fields[-1].strip(',').split(',')[0])) #UTR end(0-base): Trx start site + UTR end position in trx(0)
                UTR_end = str(int(UTR_start) + int(fields[-2].strip(',').split(',')[0])) #UTR start: UTR start + first exon length in genome
                this_UTR = fields[0] + UTR_end + curr_strand #Chrom - UTR start site - stand(+/-)
                UTR_end_new = UTR_start
            else: #No strand information
                continue
            
            if not this_UTR in scanned_3UTR_list: ##TODO: 終止コドンが同じ場合、3'UTRが最も長いものを選択すべき(Bedファイル内ですでにそのように整列してある？)
                num_exons = fields[9] #Number of exon

                #Remove shortRNA(e.g. miRNA, snoRNA...) information
                if int(num_exons) == 1:
                    if abs(int(UTR_end) - int(UTR_start)) < 200: #Remove <200 bp transripts
                        continue

                #Save 3UTR information
                write_line = [fields[0], UTR_start, UTR_end, UTR_id_new, num_exons, curr_strand] #6bed format
                #print("\t".join(write_line), end="\n", file=output_write) #print out 3UTR information
                scanned_3UTR_list[this_UTR] = write_line #Check this UTR
                raw_utr_dict[UTR_id_new] = write_line #Reserve 3UTR information as dictionary
                #num_saved += 1 #Count the number of transripts passed criteria

            else: #Already existed this 3UTR information
                UTR_end_old = scanned_3UTR_list[this_UTR][2]
                if curr_strand == '+' and UTR_end_old < UTR_end_new:
                    write_line = [fields[0], UTR_start, UTR_end, UTR_id_new, num_exons, curr_strand] #6bed format
                    scanned_3UTR_list[this_UTR] = write_line #Check this UTR
                    raw_utr_dict[UTR_id_new] = write_line #Reserve 3UTR information as dictionary
                elif curr_strand == '-' and UTR_end_old > UTR_end_new:
                    write_line = [fields[0], UTR_start, UTR_end, UTR_id_new, num_exons, curr_strand] #6bed format
                    scanned_3UTR_list[this_UTR] = write_line #Check this UTR
                    raw_utr_dict[UTR_id_new] = write_line #Reserve 3UTR information as dictionary

    #Print out 3UTR information
    for x in raw_utr_dict.keys():
        write_line = raw_utr_dict[x]
        print("\t".join(write_line), end="\n", file=output_write) #print out 3UTR information
        num_saved += 1 #Count the number of transripts passed criteria

    output_write.close()
    return raw_utr_dict
    print("Total extracted 3'UTR: " + str(num_saved))
