#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
# import pandas as pd
from Bio import SeqIO
from optparse import OptionParser


def chunkstring(string, length):
    positionlist = {}
    m = 0
    for i in range(0, len(string), int(length)):
        positionlist[m] = i + int(length)
        m += 1
    return positionlist


def splitgenome(in_genome, length):
    genome = open(in_genome, 'r')
    genomeseqs = SeqIO.to_dict(SeqIO.parse(genome, 'fasta'))
    bins = {}
    for chr in genomeseqs:
        bins[chr] = chunkstring(genomeseqs[chr], length)
    return bins


def dictxls(in_xls):
    xls = open(in_xls, 'r')
    lines = xls.readlines()
    peaks_list = []
    for line in lines:
        if line[0] != '#' and line != '\n' and 'start' not in line:
            peaks_list.append(line.rstrip('\n'))
    peak_chr = {}
    peak_start = {}
    peak_end = {}
    peak_name = {}
    m = 0
    for peak in peaks_list:
        chr = peak.split('\t')[0]
        start = peak.split('\t')[1]
        end = peak.split('\t')[2]
        name = peak.split('\t')[9]
        peak_chr[m] = chr
        peak_start[m] = start
        peak_end[m] = end
        peak_name[m] = name
        m += 1
    return peak_chr, peak_start, peak_end, peak_name


def findcommon():
    parser = OptionParser()
    parser.add_option("--genome", help="genome file")
    parser.add_option("--length", default="800", help="length of bin")
    parser.add_option("--peakfile_dirpath", help="dirpath of peakfile")
    options, args = parser.parse_args()
    bins = splitgenome(options.genome, options.length)
    binlist = []
    # binlog = open('./' + (options.peakfile_dirpath.split('/')[-1]).split('.')[0] + '.log', 'w')
    filename = options.peakfile_dirpath.split('/')[-1]
    prefix = os.path.splitext(filename)[0]
    bintype_outfile = open(prefix + '.bintype.xls', 'w')
    peak_chr, peak_start, peak_end, peak_name = dictxls(options.peakfile_dirpath)
    types = []
    bintype = {}
    for chr in bins:
        bin = bins[chr]
        for b in bin:
            type = 0
            for key in peak_name:
                if peak_chr[key] == chr:
                    if (int(bin[b]) - int(options.length)) <= int(peak_start[key]) and \
                            int(peak_start[key]) <= int(bin[b]):
                        type = 1
                        break
                    elif (int(bin[b]) - int(options.length)) <= int(peak_end[key]) and \
                            int(peak_end[key]) <= int(bin[b]):
                        type = 1
                        break
                    elif (int(bin[b]) - int(options.length)) >= int(peak_start[key]) and \
                            int(bin[b]) <= int(peak_end[key]):
                        type = 1
                        break
                    else:
                        type = 0
            types.append(type)
            binname = str(chr) + '_bin_' + str(bin[b] - int(options.length) + 1) + '_' + str(bin[b])
            binlist.append(binname)
            bintype[binname] = type
            # binlog.write(binname + ' checked' + '\n')
    for binname in bintype:
        bintype_outfile.write(binname + '\t' + str(bintype[binname]) + '\n')
    # binlog.close()
    bintype_outfile.close()


findcommon()