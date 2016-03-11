#!/usr/bin/env python3

import re

def Extract_gene_symbol_map_kfXref_file(gene_gtf_file):
    ref_dict = {}
    for line in open(gene_gtf_file, 'r'):
        if re.match("^#", line):
            continue
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
        refseq_trx_id = fields[0]
        gene_symbol = fields[1]
        refseq_trx_gene_symbol_dict[refseq_trx_id] = gene_symbol

    #for each transcript
    scanned_3UTR_list = {}
    num_saved = 0
    for line in open(gene_bed_file, 'r'):
        fields = line.rstrip().split("\t")
        chrom = fields[0]
        refseq_id = fields[3] #refseq_id
        if '_' not in fields[0]: #chr1, chr2...
            if not refseq_id in refseq_trx_gene_symbol_dict:
                gene_symbol = "NA"
                print("WARNINGS: " + gene_symbol + " does not exist...")
            else:
                gene_symbol = refseq_trx_gene_symbol_dict[refseq_id]

            UTR_id = [refseq_id, gene_symbol, fields[0], fields[5]] #refseq_id - gene_symbol - chrom - strand
            UTR_id = list(map(str, UTR_id))
            UTR_id_new = '|'.join(UTR_id)

            CDS_start = int(fields[6])
            CDS_end = int(fields[7])
            if CDS_start == CDS_end: #ncRNA (not mRNA)
                continue

            trx_start = int(fields[1])
            trx_end = int(fields[2])
            curr_strand = fields[5]
            if curr_strand == '+':
                if (trx_end - CDS_end) == 0: #3'UTR does not exist
                    continue
            elif curr_strand == '-': #3'UTR does not exist
                if (CDS_start - trx_start) == 0:
                    continue

            UTR_start = ''
            UTR_end = ''
            this_UTR = ''
            UTR_end_new = ''
            if curr_strand == '+': #Strand: +
                UTR_start = CDS_end #UTR_start (0-based)
                #this_UTR = '|'.join(map(str, [chrom, UTR_start, curr_strand]))
            elif curr_strand == '-':
                UTR_start = CDS_start #UTR_start (0-based)
                #this_UTR = '|'.join(map(str, [chrom, UTR_start, curr_strand]))
            else: #No strand information
                continue

            #if not this_UTR in scanned_3UTR_list: ##TODO: 終止コドンが同じ場合、3'UTRが最も長いものを選択すべき(Bedファイル内ですでにそのように整列してある？)
            if 0 == 0:
                num_exons = fields[9] #Number of exon

                #Remove shortRNA(e.g. miRNA, snoRNA...) information
                if int(num_exons) == 1:
                    if abs(int(trx_start) - int(trx_end)) < 200: #Remove <200 bp transripts
                        continue

                #Save 3UTR information
                blockSizes = fields[10].strip(',').split(',')
                blockStarts = fields[11].strip(',').split(',')

                write_line = []
                if curr_strand == '+': #Index => '-'
                    #Define exons with 3'UTR (index)
                    chrom_exon_st_list = [trx_start + int(blockStarts[x]) + int(blockSizes[x]) for x in range(len(blockStarts))] #exon end
                    check_block_index = len([x for x in chrom_exon_st_list if (UTR_start - x) <= 0])

                    #Extract defined exons
                    blockStarts_UTR = blockStarts[-check_block_index:]
                    blockStarts_UTR_mod = [int(x) - int(blockStarts_UTR[0]) for x in blockStarts_UTR]
                    blockSizes_UTR = blockSizes[-check_block_index:]
                    st_UTR = trx_start + int(blockStarts_UTR[0])

                    blockSizes_UTR = list(map(str, blockSizes_UTR))
                    blockStarts_UTR_mod = list(map(str, blockStarts_UTR_mod))                    

                    #Chrom - TRXstart - TRXend - name - strand - exonNumber - exonSizes - exonStarts
                    symbol = refseq_trx_gene_symbol_dict[refseq_id]
                    write_line = [chrom, st_UTR, trx_end, refseq_id, symbol, curr_strand, check_block_index, ','.join(blockSizes_UTR), ','.join(blockStarts_UTR_mod)]

                    #if check_block_index == 1:
                          

                elif curr_strand == '-': #Index => '+'
                    #Define exons with 3'UTR (index)
                    chrom_exon_st_list = [trx_start + int(blockStarts[x]) for x in range(len(blockStarts))] #exon start
                    check_block_index = len([x for x in chrom_exon_st_list if (UTR_start - x) >= 0])

                    #Extract defined exons
                    blockStarts_UTR = blockStarts[:check_block_index]
                    blockSizes_UTR = blockSizes[:check_block_index]
                    ed_UTR = trx_start + int(blockStarts_UTR[-1]) + int(blockSizes_UTR[-1])

                    blockSizes_UTR = list(map(str, blockSizes_UTR))
                    blockStarts_UTR = list(map(str, blockStarts_UTR))

                    #Chrom - TRXstart - TRXend - name - strand - exonNumber - exonSizes - exonStarts
                    symbol = refseq_trx_gene_symbol_dict[refseq_id]
                    write_line = [chrom, trx_start, ed_UTR, refseq_id, symbol, curr_strand, check_block_index, ','.join(blockSizes_UTR), ','.join(blockStarts_UTR)]

                write_line = list(map(str, write_line))
                print("\t".join(write_line), end="\n", file=output_write) #print out 3UTR information



    #            write_line = [fields[0], UTR_start, UTR_end, UTR_id_new, num_exons, curr_strand] #6bed format
    #            #print("\t".join(write_line), end="\n", file=output_write) #print out 3UTR information
    #            scanned_3UTR_list[this_UTR] = write_line #Check this UTR
    #            raw_utr_dict[UTR_id_new] = write_line #Reserve 3UTR information as dictionary
    #            #num_saved += 1 #Count the number of transripts passed criteria

    #        else: #Already existed this 3UTR information
    #            UTR_end_old = scanned_3UTR_list[this_UTR][2]
    #            if curr_strand == '+' and UTR_end_old < UTR_end_new:
    #                write_line = [fields[0], UTR_start, UTR_end, UTR_id_new, num_exons, curr_strand] #6bed format
    #                scanned_3UTR_list[this_UTR] = write_line #Check this UTR
    #                raw_utr_dict[UTR_id_new] = write_line #Reserve 3UTR information as dictionary
    #            elif curr_strand == '-' and UTR_end_old > UTR_end_new:
    #                write_line = [fields[0], UTR_start, UTR_end, UTR_id_new, num_exons, curr_strand] #6bed format
    #                scanned_3UTR_list[this_UTR] = write_line #Check this UTR
    #                raw_utr_dict[UTR_id_new] = write_line #Reserve 3UTR information as dictionary

    ##Print out 3UTR information
    #for x in raw_utr_dict.keys():
    #    write_line = raw_utr_dict[x]
    #    print("\t".join(write_line), end="\n", file=output_write) #print out 3UTR information
    #    num_saved += 1 #Count the number of transripts passed criteria

    output_write.close()
    return "NA"
    print("Total extracted 3'UTR: " + str(num_saved))
