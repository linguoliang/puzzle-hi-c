#!/usr/bin/env python
# coding: utf-8

import argparse
from array import array
import glob
import os
import pickle
import shutil
import subprocess
import sys
from multiprocessing import Pool, cpu_count
from scipy.signal import find_peaks
from scipy.ndimage import uniform_filter1d
# import matplotlib.pyplot as plt

# import argparse
import h5py
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
# In[1]:
import pandas as pd
# import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utils.alignment_input import FORMATS, normalize_contacts
import utils.convert_data as converscript
import utils.generate_fasta as gf
import utils.PuzzleHiC2JBAT as JBAT
# import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--clusters', required=True, type=int, help='Chromosomes number.')
parser.add_argument('-m', '--matrix', required=True, type=str, help='Read-level contacts: Juicer, BAM/SAM, .pairs or HiC-Pro allValidPairs.')
parser.add_argument('--input-format', choices=FORMATS, default='auto',
                    help='Input format (default: infer from filename; otherwise Juicer).')
parser.add_argument('--min-mapq', type=int, default=0,
                    help='Minimum MAPQ for both ends (default: 0; requires MAPQ fields).')
export_options = parser.add_mutually_exclusive_group(required=True)
export_options.add_argument("-j", "--juicer_tools", type=str,
                            help="juicer_tools path, required for .hic export.")
export_options.add_argument('--skip-hic', action='store_true',
                            help='Skip .hic export; no Juicer installation is needed.')
parser.add_argument('-f', '--fasta', required=True, type=str, help='Scaffold fasta file.')
parser.add_argument("-p", '--prefix', default="sample", type=str, help='Output prefix! Default: sample.')
parser.add_argument('-s', '--binsize', default=10000, type=int, help='The bin size. Default: 10000.')
parser.add_argument('-t', '--cutoff', default=0.3, type=float, help='Score cutoff, 0.25-0.5 recommended. default: 0.3.')
parser.add_argument('-i', '--init_trianglesize', default=3, type=int, help='Initial triangle size. Default: 3.')
parser.add_argument('-n', '--ncpus', default=1, type=int, help='Number of threads. Default: 1.')
parser.add_argument("-e", "--error_correction",action="store_true",help="For error correction! Default: False.")
parser.add_argument("-g", "--gap",default=100, type=int,help="The size of gap between scaffolds. Default: 100.")
# parser.add_argument("-r", "--filter",action="store_true",help="Filter! Default: False")

## 定义每个连接方式序列对应的头部和尾部contig的方向，0代表头部，1代表尾部,0和1分别代表apg文件中得正向和反向
head_dict = {0: 1, 1: 0, 2: 1, 3: 0}
end_dict = {0: 0, 1: 0, 2: 1, 3: 1}
connection_tags = ['es','ss', 'ee','se']
## 0和1分别代表apg文件中得正向和反向
conection_dict={'00':'es','01':'ee','10':'ss','11':'se'}
## 反向
AGP_HEADER=["Chromosome", "Start", "End", "Order", "Tag", "Contig_ID", "Contig_start",
                                         "Contig_end", "Orientation"]
PICKLE_PROTOCOL = pickle.HIGHEST_PROTOCOL
WRITE_BUFFER_LIMIT = 10000
GLOBAL_REPEAT_DENSITY_CONTEXT = None
REPEAT_DENSITY_CONTEXT = None
LINK_CONTEXT = None


def dump_pickle(obj):
    return np.void(pickle.dumps(obj, protocol=PICKLE_PROTOCOL))


def load_pickle(value):
    if isinstance(value, np.void):
        value = value.tobytes()
    return pickle.loads(value)


def remove_path(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)


def cleanup_paths(pattern):
    for path in glob.glob(pattern):
        remove_path(path)


def merge_files(pattern, outputfile):
    with open(outputfile, "wb") as outfile:
        for path in sorted(glob.glob(pattern)):
            with open(path, "rb") as infile:
                shutil.copyfileobj(infile, outfile, length=1024 * 1024)


def merge_tmp_re_files(prefix, outputfile, clean_tmp=False):
    merge_files(os.path.join("tmp", f"{prefix}*.re"), outputfile)
    if clean_tmp:
        cleanup_paths(os.path.join("tmp", "*"))


def generate_seq(seq, list_seg):
    seqs = []
    for seg in list_seg:
        Chr = SeqRecord(seq.seq[seg[0]:seg[1]], id="{}_{}_{}".format(seq.name, seg[0], seg[1]),
                        name="{}_{}_{}".format(seq.name, seg[0], seg[1]),
                        description="{}_{}_{}".format(seq.name, seg[0], seg[1]))
        seqs.append(Chr)
    return seqs


def split_contactmat(inputfile, outputfile, agpfile):
    split_dict = {}
    scaffold_index_dict = set(agpfile.Contig_ID)
    for item in agpfile.values:
        if item[0] != item[5]:
            if item[5] in split_dict:
                split_dict[item[5]].append([int(item[6]), int(item[7])])
            else:
                split_dict[item[5]] = [[int(item[6]), int(item[7])]]
    with open(inputfile) as HiCdata:
        with open(outputfile, 'w') as Record:
            for x in HiCdata:
                x = x.strip()
                x = x.split("\t")
                if (x[1] in scaffold_index_dict) and (x[5] in scaffold_index_dict):
                    if x[1] in split_dict:
                        pos1 = int(x[2])
                        for seg in split_dict[x[1]]:
                            if seg[0] <= pos1 <= seg[1]:
                                x[1] = "{}_{}_{}".format(x[1], seg[0], seg[1])
                                x[2] = str(pos1 - seg[0] + 1)
                                break
                    if x[5] in split_dict:
                        pos1 = int(x[6])
                        for seg in split_dict[x[5]]:
                            if seg[0] <= pos1 <= seg[1]:
                                x[5] = "{}_{}_{}".format(x[5], seg[0], seg[1])
                                x[6] = str(pos1 - seg[0] + 1)
                                break
                    Record.write("\t".join(x) + '\n')


### 读取全局的重复序列密度
def read_gloable_repeat_density(filename):
    if GLOBAL_REPEAT_DENSITY_CONTEXT is None:
        with h5py.File("tmp/repeat_dict.h5", "r") as repeat_h5:
            size_dict = load_pickle(repeat_h5["size_dict"][()])
            binsize = repeat_h5["binsize"][()]
    else:
        size_dict = GLOBAL_REPEAT_DENSITY_CONTEXT["size_dict"]
        binsize = GLOBAL_REPEAT_DENSITY_CONTEXT["binsize"]

    repeat_dict = {}
    with open(filename) as inputfile:
        str1, chr1, pos1, frag1, str2, chr2, pos2, frag2 = 0, 1, 2, 3, 4, 5, 6, 7
        for rawitem in inputfile:
            item = rawitem.strip().split(maxsplit=7)
            if (item[chr1] not in size_dict) or (item[chr1] != item[chr2]):
                continue
            #             end1=Scaffolds_len_dict[item[chr1]]
            #             end2=Scaffolds_len_dict[item[chr2]]
            #             end_pos1=end1-int(item[pos1])+1
            #             end_pos2=end2-int(item[pos2])+1
            bin1 = int(item[pos1]) // binsize
            bin2 = int(item[pos2]) // binsize
            #             end_bin1=end_pos1//binsize
            #             end_bin2=end_pos2//binsize
            if bin1 == bin2:
                if bin1 < size_dict[item[1]]:
                    if item[1] not in repeat_dict:
                        repeat_dict[item[1]] = np.zeros(size_dict[item[1]], dtype=np.int32)
                    repeat_dict[item[1]][bin1] += 1
    #                 if end_bin1<size_dict[item[1]]:
    #                     repeat_dict[item[1]]["end"][end_bin1]+=1
    with h5py.File("{}.h5".format(filename), "w") as repeat_h5:
        repeat_h5["repeat_dict"] = dump_pickle(repeat_dict)


### 简化内存使用量
###具体方法，对concat数据进行排序，获得
def return_dict_matrix(maxlength):
    # matrix_dict = {}
    # a = scaffold
    # a_len = Scaffolds_len_dict[scaffold]
    a_bin = maxlength
    #     size_dict[a]=a_bin
    #     matrix_dict={}
    #     for y in range(ord_scaffold_dict[scaffold]+1,len(Contig_ID)):
    #     print(a_len)
    # b = scaffold2
    # b_len = Scaffolds_len_dict[b]
    b_bin = maxlength
    #         size_dict[b]=b_bin
    matrix_dict = {"left_up": np.zeros((a_bin, b_bin), dtype=np.int32), "left_down": np.zeros(
        (a_bin, b_bin), dtype=np.int32), "right_up": np.zeros((a_bin, b_bin), dtype=np.int32),
                   "right_down": np.zeros((a_bin, b_bin), dtype=np.int32)}
    return matrix_dict


def reset_dict_matrix(matrix_dict):
    for matrix in matrix_dict.values():
        matrix.fill(0)


def weighted_triangle_sum(counts, row_weight, col_weight, row_mask, col_mask, trianglesize):
    row_idx = np.flatnonzero(row_mask)
    col_idx = np.flatnonzero(col_mask)
    if row_idx.size == 0 or col_idx.size == 0:
        return 0.0
    weighted = counts[np.ix_(row_idx, col_idx)].astype(np.float64, copy=False)
    weighted *= row_weight[row_idx, None]
    weighted *= col_weight[col_idx][None, :]
    rows, cols = weighted.shape
    diagonal_offset = trianglesize - rows
    total = 0.0
    for row in range(rows):
        last_col = min(cols, row + diagonal_offset + 1)
        if last_col > 0:
            total += weighted[rows - row - 1, :last_col].sum(dtype=np.float64)
    return total


def score_matrix_pair(matrix_dict, norm_weight, repeat_offset_dict, x, y, trianglesize):
    lus = weighted_triangle_sum(matrix_dict["left_up"], norm_weight[x]["start"], norm_weight[y]["start"],
                                repeat_offset_dict[x]["start"], repeat_offset_dict[y]["start"], trianglesize)
    lds = weighted_triangle_sum(matrix_dict["left_down"], norm_weight[x]["end"], norm_weight[y]["start"],
                                repeat_offset_dict[x]["end"], repeat_offset_dict[y]["start"], trianglesize)
    rus = weighted_triangle_sum(matrix_dict["right_up"], norm_weight[x]["start"], norm_weight[y]["end"],
                                repeat_offset_dict[x]["start"], repeat_offset_dict[y]["end"], trianglesize)
    rds = weighted_triangle_sum(matrix_dict["right_down"], norm_weight[x]["end"], norm_weight[y]["end"],
                                repeat_offset_dict[x]["end"], repeat_offset_dict[y]["end"], trianglesize)
    return lus, lds, rus, rds


def stable_top_k_indices(row, k):
    k = min(k, row.size)
    if k <= 0:
        return np.array([], dtype=np.intp)
    threshold = np.partition(row, row.size - k)[row.size - k]
    greater = np.flatnonzero(row > threshold)
    need = k - greater.size
    if need <= 0:
        return greater
    equal = np.flatnonzero(row == threshold)
    return np.concatenate((greater, equal[-need:]))


def top_k_sum(row, k):
    k = min(k, row.size)
    if k <= 0:
        return 0
    return np.partition(row, row.size - k)[row.size - k:].sum()


def read_raw_data(filename):
    with h5py.File("tmp/{}.h5".format("rawtemp"), "r") as rawdata:
        ord_scaffold_dict = load_pickle(rawdata["ord_scaffold_dict"][()])
    with open(filename) as inputfile:
        with open(filename + ".re", 'w') as outputfile:
            tmp_write = []
            count = 0
            #             inputfile.readline()
            for item in inputfile:
                itemlist = item.strip().split('\t', 7)
                if len(itemlist) != 8:
                    continue
                if (itemlist[1] in ord_scaffold_dict) and (itemlist[5] in ord_scaffold_dict):
                    if ord_scaffold_dict[itemlist[1]] > ord_scaffold_dict[itemlist[5]]:
                        itemlist[0:4], itemlist[4:8] = itemlist[4:8], itemlist[0:4]
                    tmp_write.append("\t".join(itemlist) + '\n')
                    count += 1
                    if count >= WRITE_BUFFER_LIMIT:
                        outputfile.writelines(tmp_write)
                        tmp_write = []
                        count = 0
            if len(tmp_write) > 0:
                outputfile.writelines(tmp_write)
                # outputfile.write("\t".join(itemlist) + '\n')


def read_repeat_density(filename):
    if REPEAT_DENSITY_CONTEXT is None:
        with h5py.File("tmp/repeat_dict.h5", "r") as repeat_h5:
            size_dict = load_pickle(repeat_h5["size_dict"][()])
            binsize = repeat_h5["binsize"][()]
            Scaffolds_len_dict = load_pickle(repeat_h5["Scaffolds_len_dict"][()])
    else:
        size_dict = REPEAT_DENSITY_CONTEXT["size_dict"]
        binsize = REPEAT_DENSITY_CONTEXT["binsize"]
        Scaffolds_len_dict = REPEAT_DENSITY_CONTEXT["Scaffolds_len_dict"]

    repeat_dict = {}
    with open(filename) as inputfile:
        str1, chr1, pos1, frag1, str2, chr2, pos2, frag2 = 0, 1, 2, 3, 4, 5, 6, 7
        for rawitem in inputfile:
            item = rawitem.strip().split(maxsplit=7)
            if (item[chr1] not in size_dict) or (item[chr1] != item[chr2]):
                continue
            end1 = Scaffolds_len_dict[item[chr1]]
            end2 = Scaffolds_len_dict[item[chr2]]
            num_pos1 = min(int(item[pos1]), end1)
            num_pos2 = min(int(item[pos2]), end2)
            end_pos1 = end1 - num_pos1 + 1
            end_pos2 = end2 - num_pos2 + 1
            bin1 = num_pos1 // binsize
            bin2 = num_pos1 // binsize
            end_bin1 = end_pos1 // binsize
            # end_bin2 = end_pos2 // binsize
            if bin1 == bin2:
                if (bin1 < size_dict[item[1]]) or (end_bin1 < size_dict[item[1]]):
                    if item[1] not in repeat_dict:
                        repeat_dict[item[1]] = {"start": np.zeros(size_dict[item[1]], dtype=np.int32),
                                                "end": np.zeros(size_dict[item[1]], dtype=np.int32)}
                    if bin1 < size_dict[item[1]]:
                        repeat_dict[item[1]]["start"][bin1] += 1
                    if end_bin1 < size_dict[item[1]]:
                        repeat_dict[item[1]]["end"][end_bin1] += 1
    with h5py.File("{}.h5".format(filename), "w") as repeat_h5:
        repeat_h5["repeat_dict"] = dump_pickle(repeat_dict)


def get_links(filename):
    if LINK_CONTEXT is None:
        with h5py.File("tmp/links.h5", "r") as links:
            repeat_offset_dict = load_pickle(links["repeat_offset_dict"][()])
            size_dict = load_pickle(links["size_dict"][()])
            binsize = links["binsize"][()]
            maxlength = links["maxlength"][()]
            trianglesize = links["trianglesize"][()]
            Scaffolds_len_dict = load_pickle(links["Scaffolds_len_dict"][()])
            ord_scaffold_dict = load_pickle(links["ord_scaffold_dict"][()])
            norm_weight = load_pickle(links["norm_weight"][()])
            # Contig_ID = load_pickle(links["Contig_ID"][()])
    else:
        repeat_offset_dict = LINK_CONTEXT["repeat_offset_dict"]
        size_dict = LINK_CONTEXT["size_dict"]
        binsize = LINK_CONTEXT["binsize"]
        maxlength = LINK_CONTEXT["maxlength"]
        trianglesize = LINK_CONTEXT["trianglesize"]
        Scaffolds_len_dict = LINK_CONTEXT["Scaffolds_len_dict"]
        ord_scaffold_dict = LINK_CONTEXT["ord_scaffold_dict"]
        norm_weight = LINK_CONTEXT["norm_weight"]
    with open("{}.txt".format(filename), "w") as links:
        with open(filename) as inputfile:
            link_write = []
            rawitem = inputfile.readline()
            if not rawitem:
                return
            item = rawitem.strip().split(maxsplit=7)
            scaffold = item[1]
            scaffold2 = item[5]
            matrix_dict = return_dict_matrix(maxlength)
            left_up = matrix_dict["left_up"]
            right_up = matrix_dict["right_up"]
            left_down = matrix_dict["left_down"]
            right_down = matrix_dict["right_down"]
            str1, chr1, pos1, frag1, str2, chr2, pos2, frag2 = 0, 1, 2, 3, 4, 5, 6, 7
            inputfile.seek(0)
            for rawitem in inputfile:
                item = rawitem.strip().split(maxsplit=7)
                if (item[chr1] not in ord_scaffold_dict) or (item[chr2] not in ord_scaffold_dict) or (
                        item[chr1] == item[chr2]):
                    continue
                tmp_scaffold = item[chr1]
                tmp_scaffold2 = item[chr2]
                # if tmp_scaffold == "868186":
                #     print(item, tmp_scaffold, tmp_scaffold2)
                if ord_scaffold_dict[tmp_scaffold] > ord_scaffold_dict[tmp_scaffold2]:
                    print("error!file not sort!")
                end1 = Scaffolds_len_dict[tmp_scaffold]
                end2 = Scaffolds_len_dict[tmp_scaffold2]
                num_pos1=min(int(item[pos1]),end1)
                num_pos2 = min(int(item[pos2]), end2)
                end_pos1 = end1 - num_pos1 + 1
                end_pos2 = end2 - num_pos2 + 1
                bin1 = num_pos1 // binsize
                bin2 = num_pos2 // binsize
                end_bin1 = end_pos1 // binsize
                end_bin2 = end_pos2 // binsize
                if (scaffold != tmp_scaffold) or (scaffold2 != tmp_scaffold2):
                    # print("Procesing")
                    x = scaffold
                    y = scaffold2
                    # if x=="868186":
                    #     print(rawitem,x,y)
                    lus, lds, rus, rds = score_matrix_pair(matrix_dict, norm_weight, repeat_offset_dict,
                                                           x, y, trianglesize)

                    # lus = (np.tril(lu[::-1, ])).sum()*(1+maxlength)*maxlength/2/lums
                    # lds = (np.tril(ld[::-1, ])[::-1, ]).sum()*(1+maxlength)*maxlength/2/ldms
                    # rus = (np.tril(ru[::-1, ])[::-1, ]).sum()*(1+maxlength)*maxlength/2/rums
                    # rds = (np.tril(rd[::-1, ])[::-1, ]).sum()*(1+maxlength)*maxlength/2/rdms
                    link_write.append("{0}\t{1}\t{2:.0f}\t{3:.0f}\t{4:.0f}\t{5:.0f}\n".format(x, y, lus, lds, rus, rds))
                    if len(link_write) >= WRITE_BUFFER_LIMIT:
                        links.writelines(link_write)
                        link_write = []
                    #                     links["{}/{}".format(x,y)]=np.array([lus,lds,rus,rds],dtype=np.int32)
                    #                     wight_matrix[ord_scaffold_dict[y],ord_scaffold_dict[x]]=wight_matrix[
                    #                         ord_scaffold_dict[x],ord_scaffold_dict[y]]=max([lus,lds,rus,rds])
                    #                     tri_weight[ord_scaffold_dict[x],ord_scaffold_dict[y]]=np.argmax([lus,lds,rus,rds])
                    scaffold = tmp_scaffold
                    scaffold2 = tmp_scaffold2
                    reset_dict_matrix(matrix_dict)

                if (tmp_scaffold == scaffold) and (tmp_scaffold2 == scaffold2):
                    size1 = size_dict[tmp_scaffold]
                    size2 = size_dict[tmp_scaffold2]
                    if bin1 < size1:
                        if bin2 < size2:
                            left_up[bin1, bin2] += 1
                        if end_bin2 < size2:
                            right_up[bin1, end_bin2] += 1
                    if end_bin1 < size1:
                        if bin2 < size2:
                            left_down[end_bin1, bin2] += 1
                        if end_bin2 < size2:
                            right_down[end_bin1, end_bin2] += 1
            x = scaffold
            y = scaffold2
            lus, lds, rus, rds = score_matrix_pair(matrix_dict, norm_weight, repeat_offset_dict,
                                                   x, y, trianglesize)
            link_write.append("{0}\t{1}\t{2:.0f}\t{3:.0f}\t{4:.0f}\t{5:.0f}\n".format(x, y, lus, lds, rus, rds))
            if link_write:
                links.writelines(link_write)


def count_links(contactdata_filename, Scaffolds_len_dict, trianglesize, clusters, average_links, binsize=10000,
                Process_num=10):
    global REPEAT_DENSITY_CONTEXT, LINK_CONTEXT
    threshold = 30
    Contig_ID = []
    maxlength = trianglesize * 10
    size_dict = {}
    for scaffold in Scaffolds_len_dict:
        if Scaffolds_len_dict[scaffold] >= binsize * trianglesize:
            Contig_ID.append(scaffold)
    if len(Contig_ID) <= clusters:
        return None, None, None, None, True
    Contig_ID.sort()
    ord_scaffold_dict = {}
    for i in range(len(Contig_ID)):
        ord_scaffold_dict[Contig_ID[i]] = i
    #     print(Contig_ID)
    #     Contig_ID_Array=Array(c.c_wchar_p,Contig_ID)
    norm_weight = {}
    # matrix_dict = {}
    repeat_dict = {}
    repeat_offset_dict = {}
    #     repeat_h5=h5py.File("/data2/luoj/hicassemble/matrix_dict.h5","w")
    wight_matrix = np.zeros((len(Contig_ID), len(Contig_ID)), dtype=np.int32)
    tri_weight = np.zeros((len(Contig_ID), len(Contig_ID)), dtype=np.int8)
    tri_weight[:] = -1
    for i in range(len(Contig_ID) - 1):
        tri_weight[i, i + 1:] = 0
    print("reading raw data!")
    ## 实现多线程
    with h5py.File("tmp/{}.h5".format("rawtemp"), "w") as rawdata:
        rawdata["ord_scaffold_dict"] = dump_pickle(ord_scaffold_dict)
    subprocess.run("split -a 3 -n l/{0} -d {1} tmp/{2};".format(Process_num,
                                                               contactdata_filename, "rawtemp"), shell=True, check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    list_temp_names = []
    for i in range(Process_num):
        list_temp_names.append("tmp/{0}{1:0>3d}".format("rawtemp", i))
    with Pool(processes=Process_num) as pool:
        pool.map(read_raw_data, list_temp_names)
    merge_tmp_re_files("rawtemp", "{0}.re".format(contactdata_filename))

    print("sort raw data!")
    subprocess.run("LC_ALL=C sort -k2,2 -k6,6 {0}.re >{0}.re.sort".format(
        contactdata_filename), shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    print("processing repeat")
    for y in range(len(Contig_ID)):
        a = Contig_ID[y]
        # a_len = Scaffolds_len_dict[a]
        a_bin = maxlength
        size_dict[a] = a_bin
        repeat_dict[a] = {"start": np.zeros(a_bin, dtype=np.int32), "end": np.zeros(a_bin, dtype=np.int32)}
        repeat_offset_dict[a] = {"start": trianglesize, "end": trianglesize}
        norm_weight[a] = {"start": trianglesize, "end": trianglesize}
    print("processing repeat files")
    ### 并行化读取
    with h5py.File("tmp/repeat_dict.h5", "w") as repeat_h5:
        repeat_h5["size_dict"] = dump_pickle(size_dict)
        repeat_h5["binsize"] = binsize
        repeat_h5["Scaffolds_len_dict"] = dump_pickle(Scaffolds_len_dict)

    REPEAT_DENSITY_CONTEXT = {
        "size_dict": size_dict,
        "binsize": binsize,
        "Scaffolds_len_dict": Scaffolds_len_dict,
    }
    try:
        with Pool(processes=Process_num) as pool:
            pool.map(read_repeat_density, list_temp_names)
    finally:
        REPEAT_DENSITY_CONTEXT = None
    for temp_name in list_temp_names:
        with h5py.File("{}.h5".format(temp_name), "r") as repeat_h5:
            repeat_dict_temp = load_pickle(repeat_h5["repeat_dict"][()])
        for x, values in repeat_dict_temp.items():
            repeat_dict[x]["start"] += values["start"]
            repeat_dict[x]["end"] += values["end"]
    cleanup_paths(os.path.join("tmp", "rawtemp*"))
    remove_path("{0}.re".format(contactdata_filename))
    print("processing repeat flags")
    for x in repeat_dict:
        for y in repeat_dict[x]:
            #             mean=repeat_dict[x][y].mean()
            ## 过滤>2倍平均数和<0.5倍平均数的数据（可能是repeat）
            ## 这里修改成30或者average_links的0.3倍
            if average_links * 0.3 > 30:
                flag = repeat_dict[x][y] >= average_links * 0.3
            else:
                flag = repeat_dict[x][y] >= 30
            # norm_weight[x][y] = average_links / (repeat_dict[x][y]+0.00000001)
            #防止0的出现
            norm_weight[x][y] = average_links / (repeat_dict[x][y] +1)
            # flag = repeat_dict[x][y] >= average_links * 0.5
            # flag2 = repeat_dict[x][y] <= average_links * 2
            repeat_offset_dict[x][y] = flag

    print("processing links")
    #     print(repeat_offset_dict[x][y])
    ###splitfile
    subprocess.run("split -a 3 -n l/{0} -d {1}.re.sort tmp/{2}".format(Process_num,
                                                                       contactdata_filename, "sortemp"), shell=True,
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with h5py.File("tmp/links.h5", "w") as links:
        links["repeat_offset_dict"] = dump_pickle(repeat_offset_dict)
        links["norm_weight"] = dump_pickle(norm_weight)
        links["size_dict"] = dump_pickle(size_dict)
        links["binsize"] = binsize
        links["trianglesize"] = trianglesize
        links["maxlength"] = maxlength
        links["Scaffolds_len_dict"] = dump_pickle(Scaffolds_len_dict)
        links["ord_scaffold_dict"] = dump_pickle(ord_scaffold_dict)
        links["Contig_ID"] = dump_pickle(Contig_ID)
    list_temp_names = []
    for i in range(Process_num):
        list_temp_names.append("tmp/{0}{1:0>3d}".format("sortemp", i))
    LINK_CONTEXT = {
        "repeat_offset_dict": repeat_offset_dict,
        "norm_weight": norm_weight,
        "size_dict": size_dict,
        "binsize": binsize,
        "trianglesize": trianglesize,
        "maxlength": maxlength,
        "Scaffolds_len_dict": Scaffolds_len_dict,
        "ord_scaffold_dict": ord_scaffold_dict,
    }
    try:
        with Pool(processes=Process_num) as pool:
            pool.map(get_links, list_temp_names)
    finally:
        LINK_CONTEXT = None
    wight_matrix_mat_sparse = {}
    for temp_name in list_temp_names:
        with open("{}.txt".format(temp_name), "r") as links:
            for item in links:
                item_list = item.strip().split("\t")
                scaffold = item_list[0]
                scaffold2 = item_list[1]
                i = ord_scaffold_dict[scaffold]
                j = ord_scaffold_dict[scaffold2]
                key = (i, j)
                data = wight_matrix_mat_sparse.get(key)
                if data is None:
                    wight_matrix_mat_sparse[key] = array("i", (int(x) for x in item_list[2:]))
                else:
                    data[0] += int(item_list[2])
                    data[1] += int(item_list[3])
                    data[2] += int(item_list[4])
                    data[3] += int(item_list[5])

    for (i, j), data in wight_matrix_mat_sparse.items():
        max_value = max(data)
        wight_matrix[i, j] = max_value
        wight_matrix[j, i] = max_value
        tri_weight[i, j] = data.index(max_value)
    for temp_name in list_temp_names:
        for tmp_path in (temp_name, "{}.txt".format(temp_name)):
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    os.remove("{0}.re.sort".format(contactdata_filename))
    del wight_matrix_mat_sparse
    wight_matrix_norm = np.zeros_like(wight_matrix, dtype=np.float16)

    #     print(wight_matrix)
    matrix_len=wight_matrix.shape[1]
    for i in range(len(wight_matrix)):
        wight_matrix[i, i] = 0
        top_count = wight_matrix[i, :].max()
        if top_count <= threshold:
            wight_matrix[i, :] = 0
        tempsum = top_k_sum(wight_matrix[i, :], min(5, matrix_len))
        if tempsum == 0:
            tempsum = 1
        wight_matrix_norm[i, :] = wight_matrix[i, :] / tempsum
    target_scaffold_len_dict = {}
    for scaffold in Contig_ID:
        target_scaffold_len_dict[scaffold] = Scaffolds_len_dict[scaffold]
    # 获取统计数据
    # with h5py.File("{0}.h5".format(contactdata_filename), "w") as matrix_h5:
    #     matrix_h5["dict"]=pickle.dumps(wight_matrix_mat_raw,protocol=PICKLE_PROTOCOL)
    #     matrix_h5["ord_scaffold_dict"]=pickle.dumps(ord_scaffold_dict,protocol=PICKLE_PROTOCOL)
    #     matrix_h5["wight_matrix"] = dump_pickle(wight_matrix)
    # np.savez("{0}.matrix".format(contactdata_filename), wight_matrix)
    return wight_matrix_norm, tri_weight, Contig_ID, target_scaffold_len_dict, False


def traverse_loop(subgraph, root):
    path = [root]
    path_dict = {root: 0}
    new_path = []
    weight = []
    while True:
        next_node = return_next_node(path_dict, list(subgraph.neighbors(path[-1])))
        weight.append(subgraph.get_edge_data(path[-1], next_node)['weight'])
        #         print(path)
        if next_node in path_dict:
            break
        else:
            path.append(next_node)
            path_dict[next_node] = 0
    index = np.argmin(weight)
    index += 1
    length = len(path)
    for i in range(length):
        new_i = (index + i) % length
        new_path.append(path[new_i])
    return new_path


# def return_cutoff(n, clusters, total=3000000000, maxcutoff=0.2):
#     sum_data = 1 / maxcutoff
#     all_pre = n / clusters
#     length = total / (clusters * 1000000)
#     for i in range(0, int(all_pre) - 1):
#         sum_data += all_pre / ((i + 1) * length)
#     return 1 / (sum_data)


def traverse(subgraph, root):
    path = [root]
    nodes = list(subgraph.neighbors(path[-1]))
    if len(nodes) != 1:
        print("error!")
        sys.exit(-1)
    path.append(nodes[0])
    nodes = list(subgraph.neighbors(path[-1]))
    while len(nodes) == 2:
        next_node = return_next_node(path, nodes)
        path.append(next_node)
        nodes = list(subgraph.neighbors(path[-1]))
    return path


def return_next_node(path, nodes):
    if nodes[0] in path:
        return nodes[1]
    else:
        return nodes[0]


def orientation(path, orientations, head_dict, end_dict):
    path_orientation = [head_dict[orientations[path[0], path[1]]]]
    for i in range(1, len(path)):
        path_orientation.append(end_dict[orientations[path[i - 1], path[i]]])
    return path_orientation


## 优化cutoff策略
def generate_grahp(score, cutoff):
    # cutoff=0.4
    # cutoff = return_cutoff(n, clusters, total)
    G = nx.Graph()
    G.add_nodes_from(range(len(score)))
    if len(score) == 2:
        G.add_edge(0, 1, weight=1)
    else:
        edges1 = set([])
        edges2 = set([])
        for i in range(len(score)):
            two_max_index = stable_top_k_indices(score[i], 2)
            for j in two_max_index:
                if i < j:
                    edges1.add((i, j))
                else:
                    edges2.add((j, i))
        union_edges = edges1 & edges2
        for edge in union_edges:
            tpscore = (score[edge[0], edge[1]] + score[edge[1], edge[0]]) / 2
            if tpscore > cutoff:
                G.add_edge(edge[0], edge[1], weight=(score[edge[0], edge[1]] + score[edge[1], edge[0]]) / 2)
    return G


def generate_agp(final_path, final_path_orientation, index, index_Scaffold_dict, Scaffold_len_Dict, iteration,
                 block_index="Block"):
    temp_list = []
    for i in range(len(final_path[index])):
        Chromosome = "{0}_{1}_{2}".format(block_index, iteration, index)
        Start = 1
        End = Scaffold_len_Dict[index_Scaffold_dict[final_path[index][i]]]
        Order = 1
        Tag = "W"
        Contig_ID = index_Scaffold_dict[final_path[index][i]]
        Contig_start = Start
        Contig_end = End
        Orientation = final_path_orientation[index][i]
        temp_data = [Chromosome, Start, End, Order, Tag, Contig_ID, Contig_start, Contig_end, Orientation]
        temp_list.append(temp_data)
    tempagp = pd.DataFrame(data=temp_list,
                           columns=AGP_HEADER)
    return tempagp

def survey_contactmat(inputfile):
    with h5py.File("tmp/convert.h5", "r") as convert:
        Scaffold_dict_list = load_pickle(convert["Scaffold_dict_list"][()])
        scaffold_index_dict = load_pickle(convert["scaffold_index_dict"][()])
        fake_chrom_dict = load_pickle(convert["fake_chrom_dict"][()])
        faker_scaffold_len_dict=load_pickle(convert["faker_scaffold_len_dict"][()])
        binsize = convert["binsize"][()]
        # convert["faker_scaffold_len_dict"] = dump_pickle(faker_scaffold_len_dict)
    # add function for correction
    correct_array={}
    # with h5py.File(f"tmp/{inputfile}.h5",'w') as stats:
    with open(inputfile) as HiCdata:
        tmp_write=[]
        count=0
        # with open(inputfile + ".re", 'w') as Record:
        for x in HiCdata:
            x = x.rstrip("\n").split("\t", 7)
            if (x[1] in scaffold_index_dict) and (x[5] in scaffold_index_dict):
                chr1index = scaffold_index_dict[x[1]]
                chr1info = Scaffold_dict_list[chr1index][x[1]]
                if str(chr1info[1]) == "0":
                    pos1 = chr1info[2] + int(x[2]) - 1
                else:
                    pos1 = chr1info[3] - int(x[2]) + 1
                chr2index = scaffold_index_dict[x[5]]
                chr2info = Scaffold_dict_list  [chr2index][x[5]]
                #                     chr2info=Scaffold_dict[x[5]]
                if str(chr2info[1]) == "0":
                    pos2 = chr2info[2] + int(x[6]) - 1
                else:
                    pos2 = chr2info[3] - int(x[6]) + 1
                chr1=fake_chrom_dict[chr1index]
                chr2=fake_chrom_dict[chr2index]
                if chr1==chr2 and faker_scaffold_len_dict[chr1]>1000000 and chr1!=x[1]:
                    distant=pos2-pos1
                    if abs(distant)>100000:
                        if chr1 not in correct_array:
                            # correct_dict[chr1] = np.zeros(faker_scaffold_len_dict[chr1] // binsize + 1,
                            #                               dtype=np.int32)
                            correct_array[chr1]=np.zeros((2,faker_scaffold_len_dict[chr1] // binsize + 1),
                                                          dtype=np.int32)
                        # correct_dict[chr1][pos1 // binsize] += 1
                        # correct_dict[chr1][pos2 // binsize] += 1
                        if abs(distant)<500000:
                            if distant<0:
                                correct_array[chr1][0][pos1 // binsize] += 1
                                correct_array[chr1][1][pos2 // binsize] += 1
                            else:
                                correct_array[chr1][1][pos1 // binsize] += 1
                                correct_array[chr1][0][pos2 // binsize] += 1
                # x[1] = chr1
                # x[5] = chr2
                # x[2] = str(pos1)
                # x[6] = str(pos2)
                # # for correction
                # tmp_write.append("\t".join(x) + '\n')
                # count += 1
                # if count>1999:
                #     Record.writelines(tmp_write)
                #     tmp_write=[]
                #     count=0
            # if len(tmp_write)>0:
            #     Record.writelines(tmp_write)
    with h5py.File(f"{inputfile}.h5",'w') as stats:
        stats["correct_array"] = dump_pickle(correct_array)


def convert_contactmat(inputfile):
    with h5py.File("tmp/convert.h5", "r") as convert:
        Scaffold_dict_list = load_pickle(convert["Scaffold_dict_list"][()])
        scaffold_index_dict = load_pickle(convert["scaffold_index_dict"][()])
        fake_chrom_dict = load_pickle(convert["fake_chrom_dict"][()])
    JBAT.convert_contact_txt(inputfile, Scaffold_dict_list, scaffold_index_dict, fake_chrom_dict)

def build_index_scaffold(Scaffold_dict):
    index_Scaffold_dict = {}
    for key in range(len(Scaffold_dict)):
        index_Scaffold_dict[key] = Scaffold_dict[key]
    return index_Scaffold_dict


def mirror_orientation_matrix(oritention):
    for i in range(len(oritention) - 1):
        target = oritention[i + 1:, i]
        target[:] = oritention[i, i + 1:]
        flip_mask = (target == 1) | (target == 2)
        if flip_mask.any():
            target[flip_mask] = 3 - target[flip_mask]
    return oritention


def buil_oritention_matrix(oritention):
    return mirror_orientation_matrix(oritention)


# def sovle_thepath(score, Scaffold_len_Dict, clusters, genome_total_size):
#     G = generate_grahp(score, len(Scaffold_len_Dict), clusters, genome_total_size)
#     allpaths = []
#     for i in nx.connected_components(G):
#         path = []
#         root = list(i)[0]
#         if len(i) > 1:
#             is_loop = True
#             for j in i:
#                 if G.degree(j) == 1:
#                     root = j
#                     path.append(j)
#                     is_loop = False
#                     break
#             if is_loop:
#                 path = traverse_loop(G.subgraph(i), root)
#             else:
#                 path = traverse(G.subgraph(i), root)
#         else:
#             is_loop = False
#             path = list(i)
#             root = path[0]
#         #         print(path)
#         allpaths.append(path)
#     final_path = []
#     for path_list in allpaths:
#         temp_error = []
#         flag_pass = False
#         if len(path_list) >= 3:
#             for i in range(1, len(path_list) - 1):
#                 if flag_pass:
#                     flag_pass = False
#                     continue
#                 if ((int(oritention[path_list[i], path_list[i - 1]]) ^ int(
#                         oritention[path_list[i], path_list[i + 1]])) & 1) == 0:
#                     temp_error.append(i)
#                     flag_pass = True
#             path_list_len = len(path_list)
#             temp_error.append(path_list_len)
#             temp_start = 0
#             for index_path_list in temp_error:
#                 final_path.append(path_list[temp_start:index_path_list])
#                 if index_path_list != path_list_len:
#                     final_path.append([path_list[index_path_list]])
#                 temp_start = index_path_list + 1
#         else:
#             final_path.append(path_list)
#     #     final_path
#     ## 将每条path的scaffold 顺序定下来
#     """
#     0 表示顺序，1表示反相互补
#     """
#     final_path_orientation = []
#     for p in final_path:
#         if len(p) == 1:
#             final_path_orientation.append([0])
#         else:
#             path_orientation = orientation(p, oritention, head_dict, end_dict)
#             final_path_orientation.append(path_orientation)
#     return final_path, final_path_orientation

def find_break_point(sig):
    osig=sig
    threshold = 0.7
    sig = sig[50:-50]
    upbounder = np.quantile(sig, q=0.995)
    lowbounder = np.quantile(sig, q=0.005)
    sig[sig > upbounder] = upbounder
    sig[sig < lowbounder] = lowbounder
    smooth_signal = uniform_filter1d(sig, size=10)
    peaks, _ = find_peaks(smooth_signal, height=np.max(smooth_signal) * threshold)
    large_peaks = []
    mini_peaks = []
    if len(peaks)>0:
        minima = np.zeros_like(peaks)
        minima_arg = np.zeros_like(peaks)
        for i, peak in enumerate(peaks):
            left = peak+50 - 40
            right = peak+50 + 1
            minima[i] = np.min(osig[left:right])
            minima_arg[i] = np.argmin(osig[left:right]) + left
        for i, peak in enumerate(peaks):
            if sig[peak] - minima[i] > 0.60:  # 调整此处的阈值来筛选大锯齿信号
                large_peaks.append(peak + 50)
                mini_peaks.append(minima_arg[i])
    if len(mini_peaks)>0:
        near_break_point=[(mini_peaks[0]+large_peaks[0])/2]
        for i in range(1,len(mini_peaks)):
            if mini_peaks[i]!=mini_peaks[i-1]:
                near_break_point.append((mini_peaks[i]+large_peaks[i])/2)
    else:
        near_break_point=[]
    return near_break_point

def generate_scaffold_info(agp_list,gap=100):
    faker_scaffold_len_dict = {}
    scaffold_index_dict = {}
    fake_chrom_dict = []
    for i in range(len(agp_list)):
        agp_list[i]["Start"] = 1
        agp_list[i]["End"] = int(agp_list[i].iloc[0, 7])
        scaffold_index_dict[agp_list[i].iloc[0, 5]] = i
        fake_chrom_dict.append(agp_list[i].iloc[0, 0])
        for j in range(1, len(agp_list[i])):
            scaffold_index_dict[agp_list[i].iloc[j, 5]] = i
            agp_list[i].iloc[j, 1] = int(agp_list[i].iloc[j - 1, 2]) + 1 + gap
            agp_list[i].iloc[j, 2] = int(agp_list[i].iloc[j, 1]) + int(agp_list[i].iloc[j, 7]) - 1
        faker_scaffold_len_dict[agp_list[i].iloc[-1, 0]] = agp_list[i].iloc[-1, 2]
    ## 构建转换矩阵
    Scaffold_dict_list = []
    for i in range(len(agp_list)):
        Scaffold_dict = {}
        for x in agp_list[i].values:
            Scaffold_dict[x[5]] = [x[7], x[8], x[1], x[2]]
        Scaffold_dict_list.append(Scaffold_dict)
    ##生成迭代的agp文件
    interation_agp = pd.concat(agp_list, ignore_index=True) if agp_list else pd.DataFrame(columns=AGP_HEADER)
    return faker_scaffold_len_dict,scaffold_index_dict,fake_chrom_dict,Scaffold_dict_list,interation_agp

def split_agp(agpfile,arr):
    templist_agp=[]
    print(arr)
    for i in range(1,len(arr)):
        rows = []
        chrom_name = agpfile.iloc[0, 0] + "_break_" + str(i)
        start = 1
        for agp_row in agpfile.iloc[arr[i-1]:arr[i], :].values:
            row = list(agp_row)
            row[0] = chrom_name
            row[1] = start
            row[2] = start + int(row[7]) - 1
            rows.append(row)
            start = row[2] + 1 + 100
        templist_agp.append(pd.DataFrame(data=rows, columns=AGP_HEADER))
    return templist_agp


def find_closest_value(arr, target):
    low = 0
    high = len(arr) - 1
    midv=0
    while low < high:
        mid = (low + high) // 2
        if arr[mid] == target:
            return mid
        elif arr[mid] > target:
            if mid - 1 >= 0 and abs(arr[mid - 1] - target) < abs(arr[mid] - target):
                closest=arr[mid - 1]
                midv= mid - 1
            else:
                return mid
            high = mid -1
        else:
            if mid + 1 <= len(arr) - 1 and abs(arr[mid + 1] - target) < abs(arr[mid] - target):
                closest=arr[mid + 1]
                midv= mid + 1
            else:
                # closest=arr[mid]
                # midv= mid
                return mid
            low = mid +1
        # 更新最接近的值
        # if closest is None or abs(arr[mid] - target) < abs(closest - target):
        #     closest = arr[mid]
        #     midv=mid
    return midv

def survey_contig(list_temp_names,Process_num,debug=False):
    with Pool(processes=Process_num) as pool:
        pool.map(survey_contactmat, list_temp_names)
    # for correction purpose
    correct_dict = {}
    for correct_dict_file_name in list_temp_names:
        # print(f"")
        with h5py.File(f"{correct_dict_file_name}.h5", 'r') as stats:
            temp_correct = load_pickle(stats["correct_array"][()])
            for chr1 in temp_correct:
                if chr1 in correct_dict:
                    correct_dict[chr1] += temp_correct[chr1]
                else:
                    correct_dict[chr1] = temp_correct[chr1]
    if debug:
        subprocess.run("mkdir -p correct_file", shell=True, check=True, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
        for tpc in correct_dict:
            sig = np.log10((correct_dict[tpc][0].T + 1) / (correct_dict[tpc][1].T + 1))
            plt.plot(sig)
            plt.savefig(f"correct_file/{tpc}_log.jpg")
            plt.cla()
        with h5py.File(f"{code}_{iteration}_correct_dict.h5", 'w') as correct_file:
            correct_file["correct_dict"] = dump_pickle(correct_dict)
    return correct_dict

def sovle_link(inputfile, outputfile, score, oritention, Scaffold_dict, Scaffold_len_Dict, iteration, agpfilename, init_agp,
               cutoff, Process_num=10, binsize=10000, error_correction=False, gap=100):
    index_Scaffold_dict = dict(enumerate(Scaffold_dict))
    mirror_orientation_matrix(oritention)
    #     length=Scaffold_len_Dict[index_Scaffold_dict[0]]
    #     for scaffold_index in Scaffold_len_Dict:
    #         length=min(length,Scaffold_len_Dict[scaffold_index])
    G = generate_grahp(score, cutoff)
    allpaths = []
    for i in nx.connected_components(G):
        path = []
        root = list(i)[0]
        if len(i) > 1:
            is_loop = True
            for j in i:
                if G.degree(j) == 1:
                    root = j
                    path.append(j)
                    is_loop = False
                    break
            if is_loop:
                path = traverse_loop(G.subgraph(i), root)
            else:
                path = traverse(G.subgraph(i), root)
        else:
            is_loop = False
            path = list(i)
            root = path[0]
        #         print(path)
        allpaths.append(path)
    final_path = []
    print(f'iteration:{iteration}')
    # connections,Chrom_Dict_for_error_corection=get_all_conections(iteration,agp_iter_name,init_agp,connections,conection_dict)
    # try:
    #     all_stath5py=h5py.File(f'{inputfile}_all_stat.h5', "w-")
    # except FileExistsError:
    #     print("File already Exist")
    # else:
    #     all_stath5py.close()
    # with h5py.File("init_contact_map.h5",'r') as h5file:
    # with h5py.File(f'{inputfile}_all_stat.h5',"r+") as h5write:
    for path_list in allpaths:
        # find_error_connection(path_list,oritention,index_Scaffold_dict,connections,Chrom_Dict_for_error_corection,h5file,h5write)
        temp_error = []
        flag_pass = False
        if len(path_list) >= 3:
            for i in range(1, len(path_list) - 1):
                if flag_pass:
                    flag_pass = False
                    continue
                if ((int(oritention[path_list[i], path_list[i - 1]]) ^ int(
                        oritention[path_list[i], path_list[i + 1]])) & 1) == 0:
                    temp_error.append(i)
                    flag_pass = True
            path_list_len = len(path_list)
            temp_error.append(path_list_len)
            temp_start = 0
            for index_path_list in temp_error:
                final_path.append(path_list[temp_start:index_path_list])
                if index_path_list != path_list_len:
                    final_path.append([path_list[index_path_list]])
                temp_start = index_path_list + 1
        else:
            final_path.append(path_list)

    #     final_path
    ## 将每条path的scaffold 顺序定下来
    """
    0 表示顺序，1表示反相互补
    """
    final_path_orientation = []
    for p in final_path:
        if len(p) == 1:
            final_path_orientation.append([0])
        else:
            path_orientation = orientation(p, oritention, head_dict, end_dict)
            final_path_orientation.append(path_orientation)
    ## 将path中的scaffold连接起来形成新的scaffold
    # agp 文件生成：Chromosome	Start	End	Order	Tag	Contig_ID	Contig_start	Contig_end	Orientation
    agp_list = []
    # block_mark = {}
    for i in range(len(final_path)):
        if len(final_path[i]) == 1:
            Chromosome = index_Scaffold_dict[final_path[i][0]]
            Start = 1
            End = Scaffold_len_Dict[index_Scaffold_dict[final_path[i][0]]]
            Order = 1
            Tag = "W"
            Contig_ID = Chromosome
            Contig_start = Start
            Contig_end = End
            Orientation = final_path_orientation[i][0]
            data = [[Chromosome, Start, End, Order, Tag, Contig_ID, Contig_start, Contig_end, Orientation]]
            tempagp = pd.DataFrame(data=data,
                                   columns=AGP_HEADER)
            agp_list.append(tempagp)
        else:
            tempagp = generate_agp(final_path, final_path_orientation, i, index_Scaffold_dict, Scaffold_len_Dict,
                                   iteration)
            agp_list.append(tempagp)
    faker_scaffold_len_dict, scaffold_index_dict, fake_chrom_dict, Scaffold_dict_list, interation_agp = generate_scaffold_info(
        agp_list,gap)
    ###生成正确的chrom
    # faker_scaffold_len_dict = {}
    # scaffold_index_dict = {}
    # fake_chrom_dict = []
    # for i in range(len(agp_list)):
    #     agp_list[i]["Start"] = 1
    #     agp_list[i]["End"] = int(agp_list[i].iloc[0, 7])
    #     scaffold_index_dict[agp_list[i].iloc[0, 5]] = i
    #     fake_chrom_dict.append(agp_list[i].iloc[0, 0])
    #     for j in range(1, len(agp_list[i])):
    #         scaffold_index_dict[agp_list[i].iloc[j, 5]] = i
    #         agp_list[i].iloc[j, 1] = int(agp_list[i].iloc[j - 1, 2]) + 1 + 100
    #         agp_list[i].iloc[j, 2] = int(agp_list[i].iloc[j, 1]) + int(agp_list[i].iloc[j, 7]) - 1
    #     faker_scaffold_len_dict[agp_list[i].iloc[-1, 0]] = agp_list[i].iloc[-1, 2]
    # ## 构建转换矩阵
    # Scaffold_dict_list = []
    # for i in range(len(agp_list)):
    #     Scaffold_dict = {}
    #     for x in agp_list[i].values:
    #         Scaffold_dict[x[5]] = [x[7], x[8], x[1], x[2]]
    #     Scaffold_dict_list.append(Scaffold_dict)
    ##生成迭代的agp文件
    # interation_agp = pd.DataFrame(data=[],
    #                               columns=["Chromosome", "Start", "End", "Order", "Tag", "Contig_ID", "Contig_start",
    #                                        "Contig_end", "Orientation"])
    # for tempagp in agp_list:
    #     interation_agp = interation_agp.append(tempagp)
    interation_agp.to_csv(agpfilename.format(iteration), sep="\t",header=False, index=False)
    subprocess.run("split -a 3 -n l/{0} -d {1} tmp/{2};".format(Process_num,
                                                               inputfile, "convertemp"), shell=True, check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if iteration>0:
        remove_path(inputfile)
    with h5py.File("tmp/convert.h5", "w") as convert:
        convert["Scaffold_dict_list"] = dump_pickle(Scaffold_dict_list)
        convert["scaffold_index_dict"] = dump_pickle(scaffold_index_dict)
        convert["fake_chrom_dict"] = dump_pickle(fake_chrom_dict)
        convert["faker_scaffold_len_dict"] = dump_pickle(faker_scaffold_len_dict)
        convert["binsize"] = binsize
    list_temp_names = []
    for i in range(Process_num):
        list_temp_names.append("tmp/{0}{1:0>3d}".format("convertemp", i))


    # with Pool(processes=Process_num) as pool:
    #     pool.map(survey_contactmat, list_temp_names)
    # # for correction purpose
    # correct_dict = {}
    # for correct_dict_file_name in list_temp_names:
    #     # print(f"")
    #     with h5py.File(f"{correct_dict_file_name}.h5", 'r') as stats:
    #         temp_correct = load_pickle(stats["correct_array"][()])
    #         for chr1 in temp_correct:
    #             if chr1 in correct_dict:
    #                 correct_dict[chr1] += temp_correct[chr1]
    #             else:
    #                 correct_dict[chr1] = temp_correct[chr1]
    # subprocess.run("mkdir -p correct_file", shell=True, check=True, stdout=subprocess.PIPE,
    #                stderr=subprocess.PIPE)
    # for tpc in correct_dict:
    #     sig = np.log10((correct_dict[tpc][0].T + 1) / (correct_dict[tpc][1].T + 1))
    #     plt.plot(sig)
    #     plt.savefig(f"correct_file/{tpc}_log.jpg")
    #     plt.cla()
    # with h5py.File(f"{code}_{iteration}_correct_dict.h5", 'w') as correct_file:
    #     correct_file["correct_dict"] = dump_pickle(correct_dict)
    if error_correction:    # break agp file
        correct_dict = survey_contig(list_temp_names, Process_num)
        tmp_agp_list=[]
        chroms = list(pd.Categorical(interation_agp.Chromosome).categories)
        for chr1 in chroms:
            if chr1 in correct_dict:
                sig=np.log10((correct_dict[chr1][0].T+1)/(correct_dict[chr1][1].T+1))
                near_break_point=find_break_point(sig)
            else:
                near_break_point=[]
            tempagp = interation_agp[interation_agp.Chromosome == chr1]
            if len(near_break_point)>0:
                arr=list(tempagp.iloc[:, 1])
                # tempagp_list=[]
                break_index=[]
                for near_point in near_break_point:
                    break_index.append(find_closest_value(arr,near_point*binsize))
                if break_index[0]!=0:
                    break_index.insert(0,0)
                break_index.append(len(arr))
                break_index=list(set(break_index))
                break_index.sort()
                tmp_agp_list+=split_agp(tempagp,break_index)
            else:
                tmp_agp_list.append(tempagp)

            # med=np.median(correct_dict[chr1])
            # plt.hlines(med, xmin=0, xmax=len(correct_dict[chr1]), colors="r")
            # plt.hlines(med * 0.2, xmin=0, xmax=len(correct_dict[chr1]), colors="g", linestyles="--")
            # plt.hlines(med * 0.1, xmin=0, xmax=len(correct_dict[chr1]), colors="g")
            # plt.plot(a.mean(),label="mean")

        faker_scaffold_len_dict, scaffold_index_dict, fake_chrom_dict, Scaffold_dict_list, interation_agp = generate_scaffold_info(
            tmp_agp_list,gap)
        interation_agp.to_csv(agpfilename.format(iteration), sep="\t",header=False, index=False)
        with h5py.File("tmp/convert.h5", "w") as convert:
            convert["Scaffold_dict_list"] = dump_pickle(Scaffold_dict_list)
            convert["scaffold_index_dict"] = dump_pickle(scaffold_index_dict)
            convert["fake_chrom_dict"] = dump_pickle(fake_chrom_dict)
            convert["faker_scaffold_len_dict"] = dump_pickle(faker_scaffold_len_dict)
            convert["binsize"] = binsize
        list_temp_names = []
        for i in range(Process_num):
            list_temp_names.append("tmp/{0}{1:0>3d}".format("convertemp", i))
    with Pool(processes=Process_num) as pool:
        pool.map(convert_contactmat, list_temp_names)
    merge_tmp_re_files("convertemp", outputfile, clean_tmp=True)
    #     convert_contactmat(inputfile,outputfile,Scaffold_dict_list,scaffold_index_dict,fake_chrom_dict)
    return faker_scaffold_len_dict


def generate_final_agp(Chrom_Dict,gap):
    Orientation2sign={0:"+",1:"-"}
    agp_list=[]
    for chrom in Chrom_Dict:
        rows = []
        start = 1
        scaffolds = Chrom_Dict[chrom]["Scaffold"]
        orientations = Chrom_Dict[chrom]["Oritention"]
        scaffold_lens = Chrom_Dict[chrom]["Scaffold_len"]
        for i in range(len(scaffolds)):
            contig_len = int(scaffold_lens[i])
            end = start + contig_len - 1
            rows.append([chrom, start, end, 2 * i + 1, "W", scaffolds[i], 1, contig_len,
                         Orientation2sign[orientations[i]]])
            if i < len(scaffolds) - 1:
                rows.append([chrom, end + 1, end + gap, 2 * (i + 1), "U", gap, "scaffold", "yes",
                             "proximity_ligation"])
                start = end + gap + 1
        if rows:
            agp_list.append([chrom, rows[-1][2], rows])
    agp_list.sort(key=lambda i:i[1],reverse=True)
    all_rows = []
    for i in range(len(agp_list)):
        chrom_name = f"scaffold_{i+1}"
        for row in agp_list[i][2]:
            row[0] = chrom_name
            all_rows.append(row)
    return pd.DataFrame(data=all_rows, columns=AGP_HEADER)



def read_init_maps(filename):
    connection_tags = ['es','ss', 'ee','se']
    with h5py.File("tmp/init_map_params.h5", "r") as init_h5:
        size_dict = load_pickle(init_h5["size_dict"][()])
        binsize = init_h5["binsize"][()]
        Scaffolds_len_dict = load_pickle(init_h5["Scaffolds_len_dict"][()])
    str1, chr1, pos1, frag1, str2, chr2, pos2, frag2 = 0, 1, 2, 3, 4, 5, 6, 7
    with h5py.File("{}_contactmap.h5".format(filename), "w") as init_array_h5:
        with open(filename) as inputfile:
            rawitem = inputfile.readline()
            if not rawitem:
                return
            item = rawitem.strip().split()
            chr1_name = item[chr1]
            chr2_name = item[chr2]
            assay_dicts={}
            for connection_tag in connection_tags:
                assay_dicts[connection_tag] = np.zeros((size_dict[chr1_name], size_dict[chr2_name]))
            inputfile.seek(0)
            for rawitem in inputfile:
                item = rawitem.strip().split()
                int_pos1=min(int(item[pos1])-1,Scaffolds_len_dict[item[chr1]])
                int_pos2=min(int(item[pos2])-1,Scaffolds_len_dict[item[chr2]])
                if (item[chr1]==chr1_name) and (item[chr2]==chr2_name):
                    # if (int_pos2// binsize ==634) and (assay_dicts['ss'].shape[1]==634):
                    #     print(filename,chr2_name,int_pos2,Scaffolds_len_dict[chr2_name])
                    #     print(rawitem)
                    assay_dicts['ss'][int_pos1//binsize,int_pos2// binsize]+=1
                    assay_dicts['se'][int_pos1 // binsize, (Scaffolds_len_dict[chr2_name]-int_pos2) // binsize] += 1
                    assay_dicts['es'][(Scaffolds_len_dict[chr1_name] - int_pos1) // binsize, int_pos2 // binsize] += 1
                    assay_dicts['ee'][(Scaffolds_len_dict[chr1_name] - int_pos1) // binsize, (Scaffolds_len_dict[chr2_name] - int_pos2) // binsize] += 1
                    if chr1_name==chr2_name:
                        if (int_pos1 // binsize!=int_pos2 // binsize):
                            assay_dicts['ss'][int_pos2 // binsize,int_pos1 // binsize] += 1
                        if (int_pos1 // binsize!= (Scaffolds_len_dict[chr2_name] - int_pos2) // binsize):
                            assay_dicts['se'][(Scaffolds_len_dict[chr2_name] - int_pos2) // binsize,int_pos1 // binsize] += 1
                        if ((Scaffolds_len_dict[chr1_name] - int_pos1) // binsize!=int_pos2 // binsize):
                            assay_dicts['es'][int_pos2 // binsize,(Scaffolds_len_dict[chr1_name] - int_pos1) // binsize] += 1
                        if ((Scaffolds_len_dict[chr1_name] - int_pos1) // binsize!=(Scaffolds_len_dict[chr2_name] - int_pos2) // binsize):
                            assay_dicts['ee'][(Scaffolds_len_dict[chr2_name] - int_pos2) // binsize,(Scaffolds_len_dict[chr1_name] - int_pos1) // binsize] += 1
                else:
                    for connection_tag in connection_tags:
                        if f'{chr1_name}/{chr2_name}/{connection_tag}' in init_array_h5:
                            print(f"read_init_maps: error in sorting contactmapps in {filename}!")
                            init_array_h5[f'{chr1_name}/{chr2_name}/{connection_tag}'][...]=assay_dicts[connection_tag]
                        else:
                            init_array_h5.create_dataset(f"{chr1_name}/{chr2_name}/{connection_tag}",data=assay_dicts[connection_tag])
                    # 初始化
                    chr1_name=item[chr1]
                    chr2_name=item[chr2]
                    assay_dicts = {}
                    for connection_tag in connection_tags:
                        assay_dicts[connection_tag] = np.zeros((size_dict[chr1_name], size_dict[chr2_name]))
                    assay_dicts['ss'][int_pos1 // binsize, int_pos2 // binsize] += 1
                    assay_dicts['se'][int_pos1 // binsize, (Scaffolds_len_dict[chr2_name] - int_pos2) // binsize] += 1
                    assay_dicts['es'][(Scaffolds_len_dict[chr1_name] - int_pos1) // binsize, int_pos2 // binsize] += 1
                    assay_dicts['ee'][(Scaffolds_len_dict[chr1_name] - int_pos1) // binsize, (
                                Scaffolds_len_dict[chr2_name] - int_pos2) // binsize] += 1
                    if chr1_name==chr2_name:
                        if (int_pos1 // binsize!=int_pos2 // binsize):
                            assay_dicts['ss'][int_pos2 // binsize,int_pos1 // binsize] += 1
                        if (int_pos1 // binsize!= (Scaffolds_len_dict[chr2_name] - int_pos2) // binsize):
                            assay_dicts['se'][(Scaffolds_len_dict[chr2_name] - int_pos2) // binsize,int_pos1 // binsize] += 1
                        if ((Scaffolds_len_dict[chr1_name] - int_pos1) // binsize!=int_pos2 // binsize):
                            assay_dicts['es'][int_pos2 // binsize,(Scaffolds_len_dict[chr1_name] - int_pos1) // binsize] += 1
                        if ((Scaffolds_len_dict[chr1_name] - int_pos1) // binsize!=(Scaffolds_len_dict[chr2_name] - int_pos2) // binsize):
                            assay_dicts['ee'][(Scaffolds_len_dict[chr2_name] - int_pos2) // binsize,(Scaffolds_len_dict[chr1_name] - int_pos1) // binsize] += 1
            for connection_tag in connection_tags:
                if f'{chr1_name}/{chr2_name}/{connection_tag}' in init_array_h5:
                    print(f"read_init_maps: error in sorting contactmapps in {filename} for the last!")
                    init_array_h5[f'{chr1_name}/{chr2_name}/{connection_tag}'][...] = assay_dicts[connection_tag]
                else:
                    init_array_h5.create_dataset(f"{chr1_name}/{chr2_name}/{connection_tag}",
                                                 data=assay_dicts[connection_tag])


    # build init contactmap
def create_init_contact_map(init_contact,Scaffolds_len_dict,Process_num,list_temp_names,size_dict_init,binsize):
    init_contact_map_h5 = h5py.File("./init_contact_map.h5", 'w')
    subprocess.run("LC_ALL=C sort -k2,2 -k6,6 {0}>{1}".format(init_contact,
                                                                    init_contact + ".sort"), shell=True, check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run("split -a 3 -n l/{0} -d {1} tmp/{2}".format(Process_num,
                                                               init_contact + ".sort", "inittemp"), shell=True, check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ### 并行化读取
    with h5py.File("tmp/init_map_params.h5", "w") as init_h5:
        init_h5["size_dict"] = dump_pickle(size_dict_init)
        init_h5["binsize"] = binsize
        init_h5["Scaffolds_len_dict"] = dump_pickle(Scaffolds_len_dict)
    with Pool(processes=Process_num) as pool:
        pool.map(read_init_maps, list_temp_names)
    for temp_name in list_temp_names:
        with h5py.File("{}_contactmap.h5".format(temp_name), "r") as init_array_h5:
            for fkey in init_array_h5.keys():
                for skey in init_array_h5[fkey].keys():
                    for tkey in init_array_h5[f'{fkey}/{skey}'].keys():
                        if f'{fkey}/{skey}/{tkey}' in init_contact_map_h5:
                            tempdata=init_contact_map_h5[f'{fkey}/{skey}/{tkey}'][...]
                            tempdata=tempdata+init_array_h5[f'{fkey}/{skey}/{tkey}'][...]
                            init_contact_map_h5[f'{fkey}/{skey}/{tkey}'][...]=tempdata
                        else:
                            init_contact_map_h5.create_dataset(f'{fkey}/{skey}/{tkey}',data=init_array_h5[f'{fkey}/{skey}/{tkey}'][...])
    os.remove(init_contact + ".sort")
    init_contact_map_h5.close()

## 生成链接路线
def generage_agp_from_scratch():
    pass

# 获取所有的连接
def get_all_conections(iteration,agp_iter_name,init_agp,connections,conection_dict):
    pd_data_list = []
    if iteration==0:
        pd_data_list.append(pd.read_csv(init_agp, names=AGP_HEADER,sep='\t',index_col=False))
    else:
        for i in range(iteration):
            temp_pd = pd.read_csv(agp_iter_name.format(i), names=AGP_HEADER,sep='\t',index_col=False)
            pd_data_list.append(temp_pd)
    Chrom_list = list(pd.Categorical(pd_data_list[-1].Chromosome).categories)
    Chrom_Dict = {}
    for i in Chrom_list:
        Chrom_Dict[i] = {}
        temp_agp = pd_data_list[-1][pd_data_list[-1].Chromosome == i]
        Chrom_Dict[i]["Scaffold"] = list(temp_agp.Contig_ID)
        Chrom_Dict[i]["Oritention"] = list(temp_agp.Orientation)
        Chrom_Dict[i]["Scaffold_len"] = list(temp_agp.Contig_end)
    for i in range(len(pd_data_list) - 2, -1, -1):
        templist = set([])
        for chrom in Chrom_Dict:
            temp_Scaffold = []
            temp_Oritention = []
            temp_Scaffold_len = []
            for j in range(len(Chrom_Dict[chrom]["Scaffold"])):
                templist.add(Chrom_Dict[chrom]["Scaffold"][j])
                temp_agp = pd_data_list[i][pd_data_list[i].Chromosome == Chrom_Dict[chrom]["Scaffold"][j]]
                #             print()
                if Chrom_Dict[chrom]["Oritention"][j] == 1:
                    temp_Scaffold.extend(list(temp_agp.Contig_ID[::-1]))
                    temp_Oritention.extend(list(1 - temp_agp.Orientation[::-1]))
                    temp_Scaffold_len.extend(list(temp_agp.Contig_end[::-1]))
                else:
                    temp_Scaffold.extend(list(temp_agp.Contig_ID))
                    temp_Oritention.extend(list(temp_agp.Orientation))
                    temp_Scaffold_len.extend(list(temp_agp.Contig_end))
            Chrom_Dict[chrom]["Scaffold"] = temp_Scaffold
            Chrom_Dict[chrom]["Oritention"] = temp_Oritention
            Chrom_Dict[chrom]["Scaffold_len"] = temp_Scaffold_len
        temchrom = set(pd_data_list[i].Chromosome)
        newset = temchrom - templist
        for k in newset:
            Chrom_Dict[k] = {}
            temp_agp = pd_data_list[i][pd_data_list[i].Chromosome == k]
            Chrom_Dict[k]["Scaffold"] = list(temp_agp.Contig_ID)
            Chrom_Dict[k]["Oritention"] = list(temp_agp.Orientation)
            Chrom_Dict[k]["Scaffold_len"] = list(temp_agp.Contig_end)
    for k in Chrom_Dict:
        for i in range(1, len(Chrom_Dict[k]["Scaffold"])):
            connection_tag = conection_dict[f'{Chrom_Dict[k]["Oritention"][i - 1]}{Chrom_Dict[k]["Oritention"][i]}']
            if f'{Chrom_Dict[k]["Scaffold"][i - 1]}/{Chrom_Dict[k]["Scaffold"][i]}/{connection_tag}' not in connections:
                connections[
                    f'{Chrom_Dict[k]["Scaffold"][i - 1]}/{Chrom_Dict[k]["Scaffold"][i]}/{connection_tag}'] = k
    return connections,Chrom_Dict

# def get_connection_contig(path_oritention,scaffold_name,connections,is_second=False):
#     if is_second:
#         path_oritention=1-path_oritention
#     if path_oritention==0:
#         contig=connections[scaffold_name]["Scaffold"][-1]
#         contig_oritention = connections[scaffold_name]["Oritention"][-1]
#     else:
#         contig = connections[scaffold_name]["Scaffold"][0]
#         contig_oritention =1- connections[scaffold_name]["Oritention"][0]
#     return contig,contig_oritention

# def find_error_connection(path,orientations,index_Scaffold_dict,Chrom_Dict,h5file,h5write):
#     # new_connection = {}
#     # 定义一下长、宽、高 width=5，length=10
#     # 最小检测范围10*binsize
#     minlength=10
#     width=5
#     if len(path)<=1:
#         return path
#     else:
#         path_oritentions = orientation(path, orientations, head_dict, end_dict)
#         for i in range(1, len(path)):
#             pre_scaffold_name=index_Scaffold_dict[path[i-1]]
#             scaffold_name=index_Scaffold_dict[path[i]]
#             if path_oritentions[i]==0:
#                 pre_contig=Chrom_Dict[pre_scaffold_name]["Scaffold"][-1]
#                 pre_contig_oritention=Chrom_Dict[pre_scaffold_name]["Oritention"][-1]
#             else:
#                 pre_contig = Chrom_Dict[pre_scaffold_name]["Scaffold"][0]
#                 pre_contig_oritention =1- Chrom_Dict[pre_scaffold_name]["Oritention"][0]
#
#             if path_oritentions[i]==0:
#                 contig=Chrom_Dict[scaffold_name]["Scaffold"][0]
#                 contig_oritention=Chrom_Dict[scaffold_name]["Oritention"][0]
#             else:
#                 contig = Chrom_Dict[scaffold_name]["Scaffold"][-1]
#                 contig_oritention =1- Chrom_Dict[scaffold_name]["Oritention"][-1]
#             connection_tag=conection_dict[f'{pre_contig_oritention}{contig_oritention}']
#             connection_tag_rev = conection_dict[f'{contig_oritention}{pre_contig_oritention}']
#             if f'{pre_contig}/{contig}/{connection_tag}' in h5file:
#                 data_for_evalue=h5file[f"{pre_contig}/{contig}/{connection_tag}"][...]
#                 data_for_pre_contig=h5file[f"{pre_contig}/{pre_contig}/{connection_tag[0]}{connection_tag[0]}"][...]
#                 data_for_contig = h5file[f"{contig}/{contig}/{connection_tag[1]}{connection_tag[1]}"][...]
#                 if min(data_for_evalue.shape)<10:
#                     break
#                 else:
#                     pre_edge=data_for_evalue[:width,:minlength]
#                     pre_contig_edge=data_for_pre_contig[:width,width:minlength+width]
#                     edge=data_for_evalue[:minlength,:width]
#                     contig_edge=data_for_contig[width:minlength+width,:width]
#             elif f'{contig}/{pre_contig}/{connection_tag_rev}' in h5file:
#                 data_for_evalue = h5file[f"{contig}/{pre_contig}/{connection_tag_rev}"][...]
#                 data_for_pre_contig = h5file[f"{pre_contig}/{pre_contig}/{connection_tag_rev[1]}{connection_tag_rev[1]}"][...]
#                 data_for_contig = h5file[f"{contig}/{contig}/{connection_tag_rev[0]}{connection_tag_rev[0]}"][...]
#                 if min(data_for_evalue.shape) < 10:
#                     break
#                 else:
#                     edge = data_for_evalue[:minlength, :width]
#                     contig_edge = data_for_pre_contig[width:minlength+width, :width]
#                     pre_edge = data_for_evalue[:width, :minlength]
#                     pre_contig_edge = data_for_contig[:width, width:minlength+width]
#             if f'{pre_contig}/{contig}/{connection_tag}' not in h5write:
#                 h5write.create_dataset(f'{pre_contig}/{contig}/{connection_tag}/pre_edge', data=pre_edge)
#                 h5write.create_dataset(f'{pre_contig}/{contig}/{connection_tag}/edge', data=edge)
#                 h5write.create_dataset(f'{pre_contig}/{contig}/{connection_tag}/pre_contig_edge', data=pre_contig_edge)
#                 h5write.create_dataset(f'{pre_contig}/{contig}/{connection_tag}/contig_edge', data=contig_edge)


def get_short_format(orig_contact):
    return normalize_contacts(orig_contact, "merged_nodups_short_format.txt", 'juicer')


# parser.add_argument("-b", "--bed", required=True, type=str, help="The bed file path!")
# parser.add_argument('-m', '--matrix', required=True, type=str, help='The matrix file path!')
if __name__ == "__main__":
    args = parser.parse_args()
    ## 判断CPU数量设定是否合理
    if args.ncpus >= cpu_count():
        print("error!cpus requrired >= ncpu!")
        sys.exit(-1)
    ###init
    ## 记录迭代次数
    iteration = 0
    ##记录聚类次数
    clusters = args.clusters
    cutoff = args.cutoff
    error_correction=args.error_correction
    # converscript = "/public/home/lgl/bin/conver_data_for_hic-Copy1.py"
    juicer_tools = args.juicer_tools
    code = args.prefix
    try:
        pair_count = normalize_contacts(args.matrix, "merged_nodups_short_format.txt",
                                        args.input_format, args.min_mapq)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f"Loaded {pair_count} contact pairs.")
    orig_contact="merged_nodups_short_format.txt"
    fastafile_name = args.fasta

    init_trianglesize = args.init_trianglesize
    growth_rate = 1.4
    trianglesize = init_trianglesize
    binsize = args.binsize
    gap=args.gap
    # check_windows = 200000
    Process_num = args.ncpus
    if not os.path.exists("./tmp"):
        os.mkdir("./tmp")
    ## for_test or for_run
    for_test = False
    # for_run = for_test

    ### init agp file
    init_agpfile_list = []
    Scaffolds_len_dict = {}

    diskcard_list = []
    split_data = {}
    # round(trianglesize*growth_rate**2)

    if for_test:
        list_chr = ["chr22", "chr15", "chr16", "chr23"]
        agp = pd.read_csv("./Carassius.auratus.build.ScafInChr.last.agp", sep="\t")
        target_agp = agp[agp.Chromosome.isin(list_chr)]
        target_agp = target_agp[target_agp.Tag == 'W']
        target_agp.Contig_end = target_agp.Contig_end.astype(int)
        for item in target_agp.values:
            Scaffolds_len_dict[item[5]] = item[7]
        contact_file = code + "_{}.txt"
        agp_iter_name = code + "_{}.agp"
        init_agp = code + "_init.agp"
        # init_contact = "{}_init.txt".format(code)
        init_contact = orig_contact
    else:
        Scaffolds_level_file = SeqIO.parse(fastafile_name, "fasta")
        for seq in Scaffolds_level_file:
            if seq.name in diskcard_list:
                pass
            elif seq.name in split_data:
                seqs = generate_seq(seq, split_data[seq.name])
                for frac_seq in seqs:
                    Scaffolds_len_dict[frac_seq.name] = len(frac_seq)
                # seqs_to_disk.extend(seqs)
            else:
                Scaffolds_len_dict[seq.name] = len(seq)
                # seqs_to_disk.append(seq)
        #     SeqIO.write(seqs_to_disk,"modified_human.fna","fasta")
        contact_file = code + "_{}.txt"
        agp_iter_name = code + "_{}.agp"
        init_agp = code + "_init.agp"
        # init_contact = "{}_init.txt".format(code)
        init_contact = orig_contact
    ## 建立初始agp 文件
    for seq in Scaffolds_len_dict:
        Chromosome = seq
        Start = 1
        End = Scaffolds_len_dict[seq]
        Order = 1
        Tag = "W"
        Contig_ID = Chromosome
        Contig_start = Start
        Contig_end = End
        Orientation = 0
        init_agpfile_list.append([Chromosome, Start, End, Order, Tag, Contig_ID, Contig_start, Contig_end, Orientation])
    init_agpfile = pd.DataFrame(data=init_agpfile_list,
                                columns=AGP_HEADER)
    init_agpfile.to_csv(init_agp, sep="\t",header=False, index=False)
    # split_contactmat(orig_contact,init_contact,init_agpfile)

    ## 分割文件，将contact matrix文件分割成多个
    subprocess.run("split -a 3 -n l/{0} -d {1} tmp/{2}".format(Process_num,
                                                               init_contact, "inittemp"), shell=True, check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    list_temp_names = []
    for i in range(Process_num):
        list_temp_names.append("tmp/{0}{1:0>3d}".format("inittemp", i))
    ##  处理repeat
    print("procesing init repeat")
    size_dict_init = {}
    repeat_dict_init = {}
    for a in Scaffolds_len_dict:
        a_len = Scaffolds_len_dict[a]
        a_bin = a_len // binsize + 1
        size_dict_init[a] = a_bin
        repeat_dict_init[a] = np.zeros(a_bin, dtype=np.int32)
    #     repeat_offset_dict[a]={"start":trianglesize,"end":trianglesize}
    print("procesing repeat files")

    ### 并行化读取
    with h5py.File("tmp/repeat_dict.h5", "w") as repeat_h5:
        # print(repeat_dict_init)
        repeat_h5["size_dict"] = dump_pickle(size_dict_init)
        repeat_h5["binsize"] = binsize
        repeat_h5["Scaffolds_len_dict"] = dump_pickle(Scaffolds_len_dict)

    GLOBAL_REPEAT_DENSITY_CONTEXT = {
        "size_dict": size_dict_init,
        "binsize": binsize,
    }
    try:
        with Pool(processes=Process_num) as pool:
            pool.map(read_gloable_repeat_density, list_temp_names)
    finally:
        GLOBAL_REPEAT_DENSITY_CONTEXT = None
    for temp_name in list_temp_names:
        with h5py.File("{}.h5".format(temp_name), "r") as repeat_h5:
            repeat_dict_temp = load_pickle(repeat_h5["repeat_dict"][()])
        for x, values in repeat_dict_temp.items():
            repeat_dict_init[x] += values
    #         repeat_dict_init[x]["end"]+=repeat_dict_temp[x]["end"]
    cleanup_paths(os.path.join("tmp", "inittemp*"))
    print("procesing repeat flags")
    all_cont_parts = [repeat_dict_init[x][:-1] for x in repeat_dict_init if repeat_dict_init[x].size > 1]
    all_cont = np.concatenate(all_cont_parts) if all_cont_parts else np.array([0], dtype=np.int32)
    # create_init_contact_map(init_contact, Scaffolds_len_dict, Process_num, list_temp_names, size_dict_init, binsize)
    # average_links = np.median(all_cont)
    ## 更改策略，仅保留70%的序列用于组装
    # reapt_flag=all_cont<average_links*1.5
    # reapt_flag1=all_cont>average_links*0.5
    # 2022-11-29
    # reapt_flag = all_cont < np.quantile(all_cont, 0.75)
    # reapt_flag1 = all_cont > np.quantile(all_cont, 0.25)
    # average_links = np.median(all_cont[reapt_flag & reapt_flag1])

    ## 采用中位数来修正
    average_links = np.median(all_cont)
    # sum_repeat=0
    # count=0
    # for x in repeat_dict_init:
    #     sum_repeat+=repeat_dict_init[x][:-1].sum()
    #     count+=len(repeat_dict_init[x][:-1])
    ## 过滤>2倍平均数和<0.5倍平均数的数据（可能是repeat）
    # average_links=sum_repeat/count
    print(f"average_links: {average_links}")
    # 减少内存
    repeat_dict_init = 0
    size_dict_init = 0





    score, oritention, Scaffold_dict, Scaffold_len_Dict, flag = count_links(init_contact, Scaffolds_len_dict,
                                                                            trianglesize,
                                                                            clusters, average_links, binsize,
                                                                            Process_num=Process_num)

    if flag:
        print("Already best assemble!")
        sys.exit(0)

    Scaffold_len_Dict = sovle_link(init_contact, contact_file.format(iteration), score, oritention, Scaffold_dict,
                                   Scaffold_len_Dict, iteration, agp_iter_name, init_agp, cutoff,
                                   Process_num=Process_num, binsize=binsize, error_correction=error_correction, gap=gap)

    for_output_dict = Scaffold_len_Dict
    iteration += 1
    while True:
        print("Processing Iter {}!".format(iteration))
        trianglesize = round(trianglesize * growth_rate)
        score, oritention, Scaffold_dict, Scaffold_len_Dict, flag = count_links(contact_file.format(iteration - 1),
                                                                                Scaffold_len_Dict, trianglesize,
                                                                                clusters,
                                                                                average_links, binsize,
                                                                                Process_num=Process_num)
        if flag:
            print("Reach the best!")
            remove_path(contact_file.format(iteration - 1))
            break
        Scaffold_len_Dict = sovle_link(contact_file.format(iteration - 1), contact_file.format(iteration), score,
                                       oritention, Scaffold_dict, Scaffold_len_Dict, iteration, agp_iter_name, init_agp,
                                       cutoff, Process_num=Process_num, binsize=binsize, error_correction=error_correction, gap=gap)
        if len(for_output_dict) >= clusters > len(Scaffold_len_Dict):
            print("Reach the best with {} Iterations".format(iteration - 2))
            iteration -= 1
            break
        for_output_dict = Scaffold_len_Dict
        iteration += 1

    # In[13]:
    scaffold_list = list(for_output_dict.keys())
    scaffold_list.sort()
    with open("{}.Chrom.sizes".format(code), 'w') as outfiles:
        for scaffold in scaffold_list:
            outfiles.write("{}\t{}\n".format(scaffold, for_output_dict[scaffold]))

    # In[18]:
    pd_data_list = []
    pd_data_list.append(pd.read_csv(init_agp, names=AGP_HEADER,sep='\t',index_col=False))
    for i in range(iteration):
        temp_pd = pd.read_csv(agp_iter_name.format(i), names=AGP_HEADER,sep='\t',index_col=False)
        pd_data_list.append(temp_pd)
    pd_group_list = [{chrom: temp_agp for chrom, temp_agp in pd_data.groupby("Chromosome", sort=False)}
                     for pd_data in pd_data_list]
    Chrom_list = list(pd.Categorical(pd_data_list[-1].Chromosome).categories)
    Chrom_Dict = {}
    for i in Chrom_list:
        Chrom_Dict[i] = {}
        temp_agp = pd_group_list[-1][i]
        Chrom_Dict[i]["Scaffold"] = list(temp_agp.Contig_ID)
        Chrom_Dict[i]["Oritention"] = list(temp_agp.Orientation)
        Chrom_Dict[i]["Scaffold_len"] = list(temp_agp.Contig_end)
    for i in range(len(pd_data_list) - 2, -1, -1):
        pd_groups = pd_group_list[i]
        templist = set([])
        for chrom in Chrom_Dict:
            temp_Scaffold = []
            temp_Oritention = []
            temp_Scaffold_len = []
            for j in range(len(Chrom_Dict[chrom]["Scaffold"])):
                templist.add(Chrom_Dict[chrom]["Scaffold"][j])
                temp_agp = pd_groups[Chrom_Dict[chrom]["Scaffold"][j]]
                #             print()
                if Chrom_Dict[chrom]["Oritention"][j] == 1:
                    temp_Scaffold.extend(list(temp_agp.Contig_ID[::-1]))
                    temp_Oritention.extend(list(1 - temp_agp.Orientation[::-1]))
                    temp_Scaffold_len.extend(list(temp_agp.Contig_end[::-1]))
                else:
                    temp_Scaffold.extend(list(temp_agp.Contig_ID))
                    temp_Oritention.extend(list(temp_agp.Orientation))
                    temp_Scaffold_len.extend(list(temp_agp.Contig_end))
            Chrom_Dict[chrom]["Scaffold"] = temp_Scaffold
            Chrom_Dict[chrom]["Oritention"] = temp_Oritention
            Chrom_Dict[chrom]["Scaffold_len"] = temp_Scaffold_len
        temchrom = set(pd_groups)
        newset = temchrom - templist
        for k in newset:
            Chrom_Dict[k] = {}
            temp_agp = pd_groups[k]
            Chrom_Dict[k]["Scaffold"] = list(temp_agp.Contig_ID)
            Chrom_Dict[k]["Oritention"] = list(temp_agp.Orientation)
            Chrom_Dict[k]["Scaffold_len"] = list(temp_agp.Contig_end)
    all_agp = generate_final_agp(Chrom_Dict,gap)
    all_agp.to_csv("./{}.agp".format(code), sep='\t',header=False,index=False)
    gf.main("./{}.agp".format(code),fastafile_name,code)

    # In[27]:
    chrom_size_dict = JBAT.get_chrom_size_from_agp(all_agp)
    with open("{}.Chrom.sizes".format(code), 'w') as outfiles:
        chrom_keys=list(chrom_size_dict.keys())
        chrom_keys.sort(key=lambda x:int(x[9:]))
        for scaffold in chrom_keys:
            outfiles.write("{}\t{}\n".format(scaffold, chrom_size_dict[scaffold]))

    if not args.skip_hic:
        # In[19]:

        fake_chrom_dict, Scaffold_dict_list, scaffold_index_dict,faker_scaffold_len_dict = JBAT.get_convert_info(all_agp)

        # In[20]:
        # subprocess.run("rm tmp/*", shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run("split -a 3 -n l/{0} -d {1} tmp/{2}".format(Process_num,
                                                                   "merged_nodups_short_format.txt", "convertemp"),
                       shell=True,
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with h5py.File("tmp/convert.h5", "w") as convert:
            convert["Scaffold_dict_list"] = dump_pickle(Scaffold_dict_list)
            convert["scaffold_index_dict"] = dump_pickle(scaffold_index_dict)
            convert["fake_chrom_dict"] = dump_pickle(fake_chrom_dict)
            convert["faker_scaffold_len_dict"] = dump_pickle(faker_scaffold_len_dict)
            convert["binsize"] = binsize
        list_temp_names = []
        for i in range(Process_num):
            list_temp_names.append("tmp/{0}{1:0>3d}".format("convertemp", i))

        with Pool(processes=Process_num) as pool:
            pool.map(convert_contactmat, list_temp_names)
        merge_tmp_re_files("convertemp", "{}.txt".format(code), clean_tmp=True)



        # In[30]:
        converscript.convert_data(chrom_size_dict,"{}.txt".format(code),"{}.txt".format(code) + ".re")
        # subprocess.run("python {0} {1} {2} {3}".format(converscript, "{}.Chrom.sizes".format(code),
        #                                                contact_file.format(iteration),
        #                                                contact_file.format(iteration) + ".re"),
        #                shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run("LC_ALL=C sort -k2,2 -k6,6 {0}>{1}".format("{}.txt".format(code)+ ".re",
                                                                        "{}.txt".format(code) + ".re.sort"),
                       shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run("{0} pre {1} {2}.hic {3}".format(juicer_tools,
                                                        "{}.txt".format(code) + ".re.sort",code,
                                                        "{}.Chrom.sizes".format(code)), shell=True, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        remove_path("{}.txt".format(code))
        remove_path("{}.txt".format(code)+ ".re")
        remove_path("{}.txt".format(code) + ".re.sort")
