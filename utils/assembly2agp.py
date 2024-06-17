import subprocess
from collections import OrderedDict
import pandas as pd
import argparse
# from multiprocessing import Pool, cpu_count
# import sys
import PuzzleHiC2JBAT as JBAT
import convert_data as converscript

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--assembly', required=True, type=str, help='assembly file.')
# parser.add_argument('-m', '--matrix', required=True, type=str, help='The matrix file path.eg: merge_nodup.txt')
# parser.add_argument("-j", "--juicer_tools",required=True,type=str,help="juicer_tools path.")
parser.add_argument("-p", '--prefix', default="sample", type=str, help='Output prefix! Default: sample.')
parser.add_argument("-g", "--gap",default=100, type=int,help="The size of gap between scaffolds. Default: 100.")
AGP_HEADER=["Chromosome", "Start", "End", "Order", "Tag", "Contig_ID", "Contig_start",
                                         "Contig_end", "Orientation"]
def seq2agp(seq,scaffoldname,idx2scffold,gap=100):
    all_agp = []
    for i in seq:
        if i < 0:
            item = idx2scffold[abs(i)]
            oritention="-"
        else:
            oritention="+"
            item = idx2scffold[abs(i)]
        Chromosome = scaffoldname
        Start = 1
        End = int(item[2])
        Order = 1
        Tag = "W"
        Contig_ID = item[0]
        Contig_start = Start
        Contig_end = End
        Orientation = oritention
        temp_data = [Chromosome, Start, End, Order, Tag, Contig_ID, Contig_start, Contig_end, Orientation]
        gap_item=[Chromosome, Start, End,Order, "U", gap, "scaffold", "yes","proximity_ligation"]
        all_agp.append(temp_data)
        all_agp.append(gap_item)
    agp = pd.DataFrame(data=all_agp[:-1], columns=AGP_HEADER)
    agp.iloc[0, 3] = 1
    agp.iloc[0, 1] = int(agp.iloc[0, 6])
    agp.iloc[0, 2] = int(agp.iloc[0, 7])
    for i in range(1, len(agp)):
        if agp.iloc[i, 4] == "W":
            agp.iloc[i, 3] = agp.iloc[i - 1, 3] + 1
            agp.iloc[i, 1] = int(agp.iloc[i - 1, 2])
            agp.iloc[i, 2] = int(agp.iloc[i - 1, 2]) + int(agp.iloc[i, 7])
        else:
            agp.iloc[i, 3] = agp.iloc[i - 1, 3] + 1
            agp.iloc[i, 1] = int(agp.iloc[i - 1, 2])
            agp.iloc[i, 2] = int(agp.iloc[i - 1, 2]) + int(agp.iloc[i, 5])
    return agp
        # temp_list.append(temp_data)
def assembly2agp(assembly,agp):
    idx2scffold = {}
    count=0
    agp_list=[]
    all_agp = pd.DataFrame(data=[], columns=AGP_HEADER)
    with open(assembly,"r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                items=line[1:].split(" ")
                idx2scffold[int(items[1])] = items
            else:
                items=line.split(" ")
                seq = [int(x) for x in items]
                count+=1
                tmp_agp=seq2agp(seq,f"scaffold_{count}",idx2scffold,gap=100)
                agp_list.append([tmp_agp.iloc[-1, 2],tmp_agp])
                # agp_contigs=pd.concat([agp_contigs,tmp_agp])
        agp_list.sort(key=lambda i: i[0], reverse=True)
        for i in range(len(agp_list)):
            agp_list[i][1].iloc[:, 0] = f"scaffold_{i + 1}"
            all_agp = pd.concat([all_agp, agp_list[i][1]])  # for pandas2
        all_agp.to_csv(f'{agp}.agp', sep="\t", header=False,index=False)
        # return all_agp


if __name__ == '__main__':
    args = parser.parse_args()
    # fasta_file_name = args.fasta
    prefix = args.prefix
    prefix.replace(".assembly","")
    prefix.replace(".agp", "")
    # nx_list = get_nx_list(args.nx)
    # contig_len_list = get_contig_length(fasta_file_name, True)
    # nx_contig = caculate_nx(nx_list, contig_len_list)
    # print_to_screen(nx_contig, contig_len_list)
    # write_to_disk(nx_contig, contig_len_list, fasta_file_name, prefix)
    # if len(sys.argv) != 3:
    #     print("{} input.agp output.assembly".format(sys.argv[0]))
    #     sys.exit(0)
    assembly2agp(args.assembly,prefix)
    # print(args.agp)

