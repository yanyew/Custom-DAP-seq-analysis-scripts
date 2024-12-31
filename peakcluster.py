#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from optparse import OptionParser
import pandas as pd


def getpeak(inf_peakfile):
    in_peak = open(inf_peakfile, 'r')
    peak_chr = {}
    peak_start = {}
    peak_end = {}
    peak_length = {}
    peak_abs_summit = {}
    chr2peak = {}
    chr_list = []
    while True:
        line = in_peak.readline()
        if not line: break
        if line.strip().split('\t')[0] != 'chr':
            name = line.rstrip('\n').split('\t')[9]
            chr = line.rstrip('\n').split('\t')[0]
            start = line.rstrip('\n').split('\t')[1]
            end = line.rstrip('\n').split('\t')[2]
            length = line.rstrip('\n').split('\t')[3]
            abs_summit = line.rstrip('\n').split('\t')[4]
            peak_chr[name] = chr
            peak_start[name] = start
            peak_end[name] = end
            peak_length[name] = length
            peak_abs_summit[name] = abs_summit
            if chr not in chr2peak:
                chr_list.append(chr)
                chr2peak[chr] = [name]
            elif chr in chr2peak:
                chr2peak[chr] += [name]
    in_peak.close()
    return peak_chr, peak_start, peak_end, peak_length, peak_abs_summit, chr2peak, chr_list


def dictbinxls(inf_binxls):
    in_bin = open(inf_binxls, 'r')
    bin_dict = {}
    bin_chr = {}
    bin_start = {}
    bin_end = {}
    chr_list = []
    chr2bin = {}
    global bin_len
    while True:
        line = in_bin.readline()
        if not line: break
        if line.strip() != '':
            name = line.rstrip('\n').split('\t')[0]
            type = line.rstrip('\n').split('\t')[1:]
            chr = name.split('_')[0]
            start = name.split('_')[2]
            end = name.split('_')[3]
            bin_len = len(type)
            bin_dict[name] = type
            bin_chr[name] = chr
            bin_start[name] = start
            bin_end[name] = end
            if chr not in chr2bin:
                chr_list.append(chr)
                chr2bin[chr] = [name]
            elif chr in chr2bin:
                chr2bin[chr] += [name]
    in_bin.close()
    return bin_dict, bin_chr, bin_start, bin_end, chr2bin, chr_list, bin_len


def find_containing_intervals(start, end, binsize):
    start_range = start // binsize
    end_range = end // binsize
    containing_intervals = []
    if start_range == end_range:
        containing_intervals.append((start_range * binsize + 1, (start_range + 1) * binsize))
    else:
        containing_intervals.append((start_range * binsize + 1, (start_range + 1) * binsize))
        for i in range(start_range + 1, end_range):
            containing_intervals.append((i * binsize + 1, (i + 1) * binsize))
        containing_intervals.append((end_range * binsize + 1, (end_range + 1) * binsize))
    return containing_intervals


# def maxlen(list):
#     max_len = 0
#     current_len = 0
#     for item in list:
#         if item == '1':
#             current_len += 1
#         else:
#             max_len = max(max_len, current_len)
#             current_len = 0
#     max_len = max(max_len, current_len)
#     return max_len


def peakcluster():
    parser = OptionParser()
    parser.add_option("--inbinfile", help="bintype file")
    parser.add_option("--peak", help="peak region file(xls)")
    parser.add_option("--binsize", default="100", help="bin size in bintype file")
    parser.add_option("--prefix", help="prefix for output peak region file(xls)")
    options, args = parser.parse_args()
    bin_dict, bin_chr, bin_start, bin_end, chr2bin, chr_list, bin_len = dictbinxls(options.inbinfile)
    print('Dict bintype file finished')
    peak_chr, peak_start, peak_end, peak_length, peak_abs_summit, chr2peak, chr_list = getpeak(options.peak)
    print('Dict peak region finished')
    excludingpeak_outfile = open(options.prefix + '.excluding.peak.xls', 'w')
    excludingpeak_outfile.write('chr\tstart\tend\tlength\tabs_summit\tpileup\t'
                                '-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname\n')
    unionpeak_outfile = open(options.prefix + '.union.peak.xls', 'w')
    unionpeak_outfile.write('chr\tstart\tend\tlength\tabs_summit\tpileup\t'
                             '-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname\n')
    corepeak_outfile = open(options.prefix + '.core.peak.xls', 'w')
    corepeak_outfile.write('chr\tstart\tend\tlength\tabs_summit\tpileup\t'
                           '-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname\n')
    dispensablepeak_outfile = open(options.prefix + '.dispensable.peak.xls', 'w')
    dispensablepeak_outfile.write('chr\tstart\tend\tlength\tabs_summit\tpileup\t'
                                  '-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname\n')
    peak2bin = {}
    for peak in peak_chr:
        chr = peak_chr[peak]
        start = peak_start[peak]
        end = peak_end[peak]
        bins = find_containing_intervals(int(start), int(end), int(options.binsize))
        for bin in bins:
            bin_name = chr + '_bin_' + str(bin[0]) + '_' + str(bin[1])
            if peak not in peak2bin:
                peak2bin[peak] = [bin_dict[bin_name]]
            elif peak in peak2bin:
                peak2bin[peak].append(bin_dict[bin_name])
    for peak in peak2bin:
        chr = peak_chr[peak]
        start = peak_start[peak]
        end = peak_end[peak]
        length = peak_length[peak]
        abs_summit = peak_abs_summit[peak]
        df = pd.DataFrame(peak2bin[peak])
        contains_one = (df == '1').any().any()
        if contains_one:
            unionpeak_outfile.write(chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(length) +
                                     '\t' + str(abs_summit) + '\t' + '\t' + '\t' + '\t' + '\t' + peak + '\n')
            all_ones = df.eq('1').any().all()
            if all_ones:
                corepeak_outfile.write(chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(length) +
                                       '\t' + str(abs_summit) + '\t' + '\t' + '\t' + '\t' + '\t' + peak + '\n')
            else:
                dispensablepeak_outfile.write(chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(length) +
                                              '\t' + str(abs_summit) + '\t' + '\t' + '\t' + '\t' + '\t' + peak + '\n')
        else:
            excludingpeak_outfile.write(chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(length) +
                                        '\t' + str(abs_summit) + '\t' + '\t' + '\t' + '\t' + '\t' + peak + '\n')
    excludingpeak_outfile.close()
    unionpeak_outfile.close()
    corepeak_outfile.close()
    dispensablepeak_outfile.close()


peakcluster()