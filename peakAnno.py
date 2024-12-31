#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys
import os
import pandas as pd
from optparse import OptionParser


def dict_gff(infile):
    inf_gff = open(infile, 'r')
    gene_chr = {}
    gene_start = {}
    gene_end = {}
    gene_strand = {}
    exon_region = {}
    chr2gene = {}
    while True:
        line = inf_gff.readline()
        if not line: break
        if line[0] != '#' and len(line.split('\t')) == 9:
            chr = line.split('\t')[0]
            type = line.split('\t')[2]
            start = line.split('\t')[3]
            end = line.split('\t')[4]
            strand = line.split('\t')[6]
            attributes = line.split('\t')[8].rstrip('\n')
            attributes_list = attributes.split(';')
            attributes_dict = {}
            if type == 'mRNA':
                for a in attributes_list:
                    if a != '':
                        attributes_dict[a.split('=')[0]] = a.split('=')[1]
                ID = attributes_dict['ID']
                gene_chr[ID] = chr
                gene_start[ID] = start
                gene_end[ID] = end
                gene_strand[ID] = strand
                if chr not in chr2gene:
                    chr2gene[chr] = [ID]
                elif chr in chr2gene:
                    chr2gene[chr].append(ID)
            elif type == 'exon':
                for a in attributes_list:
                    if a != '':
                        attributes_dict[a.split('=')[0]] = a.split('=')[1]
                parent = attributes_dict['Parent']
                if parent not in exon_region:
                    exon_region[parent] = [[int(start), int(end)]]
                elif parent in exon_region:
                    exon_region[parent].append([int(start), int(end)])
    inf_gff.close()
    return gene_chr, gene_start, gene_end, gene_strand, exon_region, chr2gene


def dict_peak(infile):
    inf_peak = open(infile, 'r')
    peak_chr = {}
    peak_start = {}
    peak_end = {}
    peak_line = {}
    while True:
        line = inf_peak.readline()
        if not line: break
        if line[0] != '#' and len(line.split('\t')) == 10 and line.split('\t')[0] != 'chr':
            chr = line.split('\t')[0]
            start = line.split('\t')[1]
            end = line.split('\t')[2]
            name = line.rstrip('\n').split('\t')[9]
            peak_chr[name] = chr
            peak_start[name] = start
            peak_end[name] = end
            peak_line[name] = line.rstrip('\n')
    inf_peak.close()
    return peak_chr, peak_start, peak_end, peak_line


def is_overlapping(interval1, interval2):
    if interval1[1] < interval2[0] or interval2[1] < interval1[0]:
        return False
    else:
        return True


def is_interval_overlapping_with_list(interval, interval_list):
    for i in interval_list:
        if is_overlapping(interval, i):
            return True
    return False


def peak_anno():
    parser = OptionParser()
    parser.add_option("--peak", help="input peak file")
    parser.add_option("--gff", help="input genome annotation file")
    parser.add_option("--upstream", default="1000", help="length of upstream region from TSS (Promoter)")
    parser.add_option("--downstream", default="1000", help="length of downstream region from TTS")
    options, args = parser.parse_args()
    gene_chr, gene_start, gene_end, gene_strand, exon_region, chr2gene = dict_gff(options.gff)
    peak_chr, peak_start, peak_end, peak_line = dict_peak(options.peak)
    prefix = os.path.splitext(options.peak)[0]
    outf = open(prefix + '.anno.xls', 'w')
    outf.write(
        'chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname\t'
        'gene\tgene_strand\tposition\n')
    for peak in peak_chr:
        count = 0
        chr = peak_chr[peak]
        if chr in chr2gene:
            for gene in chr2gene[chr]:
                if gene_strand[gene] == '+':
                    gene_promoter = int(gene_start[gene]) - int(options.upstream)
                    gene_downstream = int(gene_end[gene]) + int(options.downstream)
                    if int(peak_end[peak]) >= gene_promoter and int(peak_end[peak]) < int(gene_start[gene]):
                        outf.write(peak_line[peak] + '\t' + gene + '\t' + gene_strand[gene] + '\t' + 'Promoter' + '\n')
                        count += 1
                    elif int(peak_start[peak]) > int(gene_end[gene]) and int(peak_start[peak]) <= gene_downstream:
                        outf.write(peak_line[peak] + '\t' + gene + '\t' + gene_strand[gene] + '\t' + 'Downstream' + '\n')
                        count += 1
                    elif int(peak_start[peak]) <= int(gene_end[gene]) and int(peak_end[peak]) >= int(gene_start[gene]):
                        overlap_with_exon = \
                            is_interval_overlapping_with_list([int(peak_start[peak]), int(peak_end[peak])], exon_region[gene])
                        if overlap_with_exon:
                            outf.write(peak_line[peak] + '\t' + gene + '\t' + gene_strand[gene] + '\t' + 'Exon' + '\n')
                            count += 1
                        else:
                            outf.write(peak_line[peak] + '\t' + gene + '\t' + gene_strand[gene] + '\t' + 'Intron' + '\n')
                            count += 1
                elif gene_strand[gene] == '-':
                    gene_promoter = int(gene_end[gene]) + int(options.upstream)
                    gene_downstream = int(gene_start[gene]) - int(options.downstream)
                    if int(peak_end[peak]) >= gene_downstream and int(peak_end[peak]) < int(gene_start[gene]):
                        outf.write(peak_line[peak] + '\t' + gene + '\t' + gene_strand[gene] + '\t' + 'Downstream' + '\n')
                        count += 1
                    elif int(peak_start[peak]) > int(gene_end[gene]) and int(peak_start[peak]) <= gene_promoter:
                        outf.write(
                            peak_line[peak] + '\t' + gene + '\t' + gene_strand[gene] + '\t' + 'Promoter' + '\n')
                        count += 1
                    elif int(peak_start[peak]) <= int(gene_end[gene]) and int(peak_end[peak]) >= int(gene_start[gene]):
                        overlap_with_exon = \
                            is_interval_overlapping_with_list([int(peak_start[peak]), int(peak_end[peak])], exon_region[gene])
                        if overlap_with_exon:
                            outf.write(peak_line[peak] + '\t' + gene + '\t' + gene_strand[gene] + '\t' + 'Exon' + '\n')
                            count += 1
                        else:
                            outf.write(peak_line[peak] + '\t' + gene + '\t' + gene_strand[gene] + '\t' + 'Intron' + '\n')
                            count += 1
        if count == 0:
            outf.write(peak_line[peak] + '\t-\t-\t' + 'Intergenic' + '\n')
    outf.close()

    df = pd.read_csv(prefix + '.anno.xls', sep="\t", header=0)
    df1 = df.drop_duplicates(subset=['name'], keep='first')
    df2 = df1.reset_index(drop=True)
    df2.to_csv(prefix + '.anno.primary.xls', sep='\t', index=False, header=True)

    count = df2['position'].value_counts()
    count_df = pd.DataFrame(count)
    count_df.to_csv(prefix + '.anno.primary.postioncount.xls', sep='\t')


peak_anno()