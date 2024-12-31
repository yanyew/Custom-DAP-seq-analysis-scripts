#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from optparse import OptionParser
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def file_to_list(file_path):
    data_list = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            data_list.append(line)
    return data_list


def correlation():
    parser = OptionParser()
    parser.add_option("--filelist", help="the path list file of files for correlation")
    parser.add_option("--out", help="the prefix of output file")
    options, args = parser.parse_args()

    file_paths = file_to_list(options.filelist)
    sets = []
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            gene_set = set(file.read().splitlines())
            sets.append(gene_set)
    correlation_matrix = pd.DataFrame(index=file_paths, columns=file_paths)
    for i in range(len(file_paths)):
        for j in range(i, len(file_paths)):
            set1 = sets[i]
            set2 = sets[j]
            intersection = len(set1 & set2)
            union = len(set1 | set2)
            correlation = intersection / union
            correlation_matrix.iloc[i, j] = correlation
            correlation_matrix.iloc[j, i] = correlation
    correlation_matrix.to_csv(options.out + '.correlation_matrix.xls', sep='\t')
    sns.heatmap(correlation_matrix, annot=True, fmt='.2f', cmap='coolwarm')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.savefig(options.out + '.heatmap.svg', format='svg')


correlation()