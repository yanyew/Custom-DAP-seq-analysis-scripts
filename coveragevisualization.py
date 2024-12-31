#!/usr/bin/env python3
# -*- coding:utf-8 -*-
from optparse import OptionParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig
import gffutils
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import os


def file_to_list(file_path):
    data_list = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            data_list.append(line)
    return data_list


def plot_transcript_structure(ax, gff_file, chromosome, start, end, fig_length):
    prefix = os.path.splitext(gff_file)[0]
    if os.path.exists(prefix + '.db'):
        db = gffutils.FeatureDB(prefix + '.db')
    else:
        db = gffutils.create_db(gff_file, prefix + '.db')

    transcripts = db.region(region=(chromosome, start, end), featuretype='mRNA')
    for transcript in transcripts:
        if transcript.strand == '+':
            arrow_rect = patches.FancyArrowPatch((transcript.start, fig_length * 0.5), (transcript.end, fig_length * 0.5),
                                                 facecolor='black', arrowstyle='simple', mutation_scale=20, edgecolor='none', linewidth= fig_length * 0.5)
            ax.add_patch(arrow_rect)
        elif transcript.strand == '-':
            arrow_rect = patches.FancyArrowPatch((transcript.end, fig_length * 0.5), (transcript.start, fig_length * 0.5),
                                                 facecolor='black', arrowstyle='simple', mutation_scale=20, edgecolor='none', linewidth= fig_length * 0.5)
            ax.add_patch(arrow_rect)
        for feature in db.children(transcript, level=1, featuretype='CDS'):
            feature_start = max(feature.start, start)
            feature_end = min(feature.end, end)
            length = feature_end - feature_start
            rect = Rectangle((feature_start, fig_length * 0.35), length, fig_length * 0.3, facecolor='black', edgecolor='none')
            ax.add_patch(rect)
        transcript_id = transcript.attributes['ID'][0]
        transcript_center = (transcript.start + transcript.end) / 2
        ax.text(transcript_center, fig_length * 0.8, transcript_id, ha='center', va='center')
    ax.axhline(0, color='lightgray', linewidth=fig_length * 0.5)
    ax.set_xlim(start, end)
    ax.axis('off')


def plot_bed(ax, bed_file, chromosome, start, end, fig_length):
    intervals = []
    with open(bed_file) as bed:
        for line in bed:
            fields = line.split('\t')
            if fields[0] == chromosome:
                interval_start = int(fields[1])
                interval_end = int(fields[2])
                if interval_end >= start and interval_start <= end:
                    intervals.append((interval_start, interval_end))
    for i, (interval_start, interval_end) in enumerate(intervals):
        rect = Rectangle((interval_start, fig_length * 0.35), interval_end - interval_start, fig_length * 0.35, facecolor='green', edgecolor='none')
        ax.add_patch(rect)
    ax.axis('off')
    ax.axhline(0, color='lightgray', linewidth=fig_length * 0.5)


def plot_bigwig(ax, bigwig_file, chromosome, start, end, limit, fig_length):
    bw = pyBigWig.open(bigwig_file)
    values = bw.values(chromosome, start, end)
    ax.bar(range(start, end), values, color='blue', align='center', edgecolor='none')
    # values = np.array(values)
    # ax.fill_between(range(start, end), values, 0, where=(values > 0), interpolate=False, step=None, color='blue', data=None)
    # ax.fill_between(range(start, end), values, 0, interpolate=False, step=None, color='blue', data=None)
    ax.axhline(0, color='lightgray', linewidth= fig_length * 0.5)
    bw.close()
    ax.axis('off')
    ax.set_ylim(0, limit)


# def plot_bedgraph(ax, bedgraph_file, chromosome, start, end, limit, fig_length):
#     values = []
#     with open(bedgraph_file) as bedgraph:
#         for line in bedgraph:
#             fields = line.split('\t')
#             if fields[0] == chromosome:
#                 interval_start = int(fields[1])
#                 interval_end = int(fields[2])
#                 value = float(fields[3])
#                 if interval_end >= start and interval_start <= end:
#                     values.extend([value] * (min(interval_end, end) - max(interval_start, start)))
#     ax.bar(range(start, end), values, color='blue', align='center', edgecolor='none')
#     # values = np.array(values)
#     # ax.fill_between(range(start, end), values, 0, where=(values > 0), interpolate=False, step=None, color='blue', data=None)
#     # ax.fill_between(range(start, end), values, 0, interpolate=False, step=None, color='blue', data=None)
#     ax.axhline(0, color='gray', linewidth= fig_length * 0.5)
#     ax.axis('off')
#     ax.set_ylim(0, limit)
def plot_bedgraph(ax, bedgraph_file, chromosome, start, end, limit, fig_length):
    with open(bedgraph_file) as bedgraph:
        for line in bedgraph:
            fields = line.split('\t')
            if fields[0] == chromosome:
                interval_start = float(fields[1])
                interval_end = float(fields[2])
                value = float(fields[3])
                if interval_end >= start and interval_start + 1 <= end:
                    rect = Rectangle((interval_start, 0), interval_end - interval_start, value, facecolor='blue',
                                     edgecolor='none')
                    ax.add_patch(rect)
    ax.axhline(0, color='lightgray', linewidth= fig_length * 0.5)
    ax.axis('off')
    ax.set_xlim(start, end)
    ax.set_ylim(0, limit)


def plot_multiple_tracks():
    parser = OptionParser()
    parser.add_option("--filetype", help="the type of coverage file, bigwig or bedgraph")
    parser.add_option("--filelist", help="the list file of coverage file")
    parser.add_option("--gff", help="gff file")
    parser.add_option("--bed", help="bed file")
    parser.add_option("--region", help="target region for plot, e.g. chr1_0_100")
    parser.add_option("--limit", default="10000", help="limit for y axis, default: 1000")
    parser.add_option("--prefix", default="output", help="prefix for output graph file, default: output")
    parser.add_option('--figsize', help='the size of subplot for each track (LengthxWidth) [default=5x1]',
                      default='5x1')
    parser.add_option("--outformat", default="svg", help="format for output graph file, default: svg")
    options, args = parser.parse_args()

    file_list = file_to_list(options.filelist)
    chromosome = (options.region).split('_')[0]
    start = int((options.region).split('_')[1])
    end = int((options.region).split('_')[2])
    limit = int(options.limit)
    fig_width = float((options.figsize).split('x')[0])
    fig_length = float((options.figsize).split('x')[1])

    global fig_height
    if options.gff != '' and options.bed != '':
        fig, axs = plt.subplots(len(file_list) + 2, 1, figsize=(fig_width, fig_length * (len(file_list) + 2)), sharex=True)
        fig.subplots_adjust(hspace=0)
        plot_transcript_structure(axs[0], options.gff, chromosome, start, end, fig_length)
        axs[0].set_xlim(start, end)

        if options.filetype == 'bigwig':
            for i, file in enumerate(file_list):
                plot_bigwig(axs[i + 1], file, chromosome, start, end, limit, fig_length)
                axs[i + 1].set_xlim(start, end)
                prefix = os.path.splitext(file)[0]
                axs[i + 1].text(start - fig_width * 0.1, fig_length, prefix, ha='right', va='center', fontsize=6)

        elif options.filetype == 'bedgraph':
            for i, file in enumerate(file_list):
                plot_bedgraph(axs[i + 1], file, chromosome, start, end, limit, fig_length)
                axs[i + 1].set_xlim(start, end)
                prefix = os.path.splitext(file)[0]
                axs[i + 1].text(start - fig_width * 0.1, fig_length, prefix, ha='right', va='center', fontsize=6)

        plot_bed(axs[len(file_list) + 1], options.bed, chromosome, start, end, fig_length)
        axs[len(file_list) + 1].set_xlim(start, end)
        fig.savefig(options.prefix + '.' + options.region + '.' + options.limit + '.' + options.outformat)

    elif options.gff == '' and options.bed != '':
        fig, axs = plt.subplots(len(file_list) + 1, 1, figsize=(fig_width, fig_length * (len(file_list) + 1)), sharex=True)
        fig.subplots_adjust(hspace=0)
        if options.filetype == 'bigwig':
            for i, file in enumerate(file_list):
                plot_bigwig(axs[i], file, chromosome, start, end, limit, fig_length)
                axs[i].set_xlim(start, end)
                prefix = os.path.splitext(file)[0]
                axs[i].text(start - fig_width * 0.1, fig_length, prefix, ha='right', va='center', fontsize=6)

        elif options.filetype == 'bedgraph':
            for i, file in enumerate(file_list):
                plot_bedgraph(axs[i], file, chromosome, start, end, limit, fig_length)
                axs[i].set_xlim(start, end)
                prefix = os.path.splitext(file)[0]
                axs[i].text(start - fig_width * 0.1, fig_length, prefix, ha='right', va='center', fontsize=6)

        plot_bed(axs[len(file_list)], options.bed, chromosome, start, end, fig_length)
        axs[len(file_list)].set_xlim(start, end)
        fig.savefig(options.prefix + '.' + options.region + '.' + options.limit + '.' + options.outformat)

    elif options.gff != '' and options.bed == '':
        fig, axs = plt.subplots(len(file_list) + 1, 1, figsize=(fig_width, fig_length * (len(file_list) + 1)), sharex=True)
        fig.subplots_adjust(hspace=0)
        plot_transcript_structure(axs[0], options.gff, chromosome, start, end, fig_length)
        axs[0].set_xlim(start, end)

        if options.filetype == 'bigwig':
            for i, file in enumerate(file_list):
                plot_bigwig(axs[i + 1], file, chromosome, start, end, limit, fig_length)
                axs[i + 1].set_xlim(start, end)
                prefix = os.path.splitext(file)[0]
                axs[i + 1].text(start - fig_width * 0.1, fig_length, prefix, ha='right', va='center', fontsize=6)

        elif options.filetype == 'bedgraph':
            for i, file in enumerate(file_list):
                plot_bedgraph(axs[i + 1], file, chromosome, start, end, limit, fig_length)
                axs[i + 1].set_xlim(start, end)
                prefix = os.path.splitext(file)[0]
                axs[i + 1].text(start - fig_width * 0.1, fig_length, prefix, ha='right', va='center', fontsize=6)

        fig.savefig(options.prefix + '.' + options.region + '.' + options.limit + '.' + options.outformat)

    elif options.gff == '' and options.bed == '':
        fig, axs = plt.subplots(len(file_list), 1, figsize=(fig_width, fig_length * len(file_list)), sharex=True)
        fig.subplots_adjust(hspace=0)

        if options.filetype == 'bigwig':
            for i, file in enumerate(file_list):
                plot_bigwig(axs[i], file, chromosome, start, end, limit, fig_length)
                axs[i].set_xlim(start, end)
                prefix = os.path.splitext(file)[0]
                axs[i].text(start - fig_width * 0.1, fig_length, prefix, ha='right', va='center', fontsize=6)

        elif options.filetype == 'bedgraph':
            for i, file in enumerate(file_list):
                plot_bedgraph(axs[i + 1], file, chromosome, start, end, limit, fig_length)
                axs[i].set_xlim(start, end)
                prefix = os.path.splitext(file)[0]
                axs[i].text(start - fig_width * 0.1, fig_length, prefix, ha='right', va='center', fontsize=6)

        fig.savefig(options.prefix + '.' + options.region + '.' + options.limit + '.' + options.outformat)


plot_multiple_tracks()