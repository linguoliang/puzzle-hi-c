import pandas as pd
import argparse
import sys

import numpy as np
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-p", '--prefix', default="Nx", type=str, help='Prefix! default:Nx')
parser.add_argument('-f', '--fasta', required=True, type=str, help='Scaffold fasta file!')
parser.add_argument('-N', '--nx', default="50,90,95,98,99", type=str, help='Nx list! default:50,90,95,98,99')

AGP_HEADER=["Chromosome", "Start", "End", "Order", "Tag", "Contig_ID", "Contig_start",
                                         "Contig_end", "Orientation"]
def get_convert_info(all_agp):
    fake_chrom_dict = list(pd.Categorical(all_agp.Chromosome).categories)
    Scaffold_dict_list = []
    scaffold_index_dict = {}
    faker_scaffold_len_dict={}
    for i in range(len(fake_chrom_dict)):
        temp_agp = all_agp[all_agp.Chromosome == fake_chrom_dict[i]]
        temp_agp=temp_agp[temp_agp.Tag == "W"]
        for sca in list(temp_agp.Contig_ID):
            scaffold_index_dict[sca] = i
        Scaffold_dict = {}
        for x in temp_agp.values:
            Scaffold_dict[x[5]] = [int(x[7]), x[8], x[1], x[2]]
        Scaffold_dict_list.append(Scaffold_dict)
        faker_scaffold_len_dict[temp_agp.iloc[-1, 0]] = temp_agp.iloc[-1, 2]
    return fake_chrom_dict, Scaffold_dict_list, scaffold_index_dict,faker_scaffold_len_dict

def convert_contact_txt(inputfile,Scaffold_dict_list,scaffold_index_dict,fake_chrom_dict):
    with open(inputfile) as HiCdata:
        tmp_write=[]
        count=0
        with open(inputfile + ".re", 'w') as Record:
            for x in HiCdata:
                x = x.strip()
                x = x.split("\t")
                if (x[1] in scaffold_index_dict) and (x[5] in scaffold_index_dict):
                    chr1index = scaffold_index_dict[x[1]]
                    chr1info = Scaffold_dict_list[chr1index][x[1]]
                    if str(chr1info[1]) == "0" or str(chr1info[1]) == "+":
                        pos1 = chr1info[2] + int(x[2]) - 1
                    else:
                        pos1 = chr1info[3] - int(x[2]) + 1
                    chr2index = scaffold_index_dict[x[5]]
                    chr2info = Scaffold_dict_list[chr2index][x[5]]
                    #                     chr2info=Scaffold_dict[x[5]]
                    if str(chr2info[1]) == "0" or str(chr2info[1]) == "+":
                        pos2 = chr2info[2] + int(x[6]) - 1
                    else:
                        pos2 = chr2info[3] - int(x[6]) + 1
                    chr1=fake_chrom_dict[chr1index]
                    chr2=fake_chrom_dict[chr2index]
                    x[1] = chr1
                    x[5] = chr2
                    x[2] = str(pos1)
                    x[6] = str(pos2)
                    # for correction
                    tmp_write.append("\t".join(x) + '\n')
                    count += 1
                    if count>1999:
                        Record.writelines(tmp_write)
                        tmp_write=[]
                        count=0
            if len(tmp_write)>0:
                Record.writelines(tmp_write)

def convert_to_supter_scaffold(agp):
    agp_contigs = agp[agp.Tag == "W"]
    agp_contigs.Chromosome="assembly"
    agp_contigs.iloc[0, 3] = 1
    agp_contigs.iloc[0,1]=agp_contigs.iloc[0,6]
    agp_contigs.iloc[0, 2] = agp_contigs.iloc[0, 7]
    for i in range(1,len(agp_contigs)):
        agp_contigs.iloc[i,3]=agp_contigs.iloc[i-1, 3]+1
        agp_contigs.iloc[i, 1] = agp_contigs.iloc[i-1, 2]
        agp_contigs.iloc[i, 2] = agp_contigs.iloc[i - 1, 2]+ agp_contigs.iloc[i, 7]
    return agp_contigs

if __name__ == "__main__":
    args = parser.parse_args()
    fasta_file_name = args.fasta
    prefix = args.prefix
    nx_list = get_nx_list(args.nx)
    contig_len_list = get_contig_length(fasta_file_name, True)
    nx_contig = caculate_nx(nx_list, contig_len_list)
    print_to_screen(nx_contig,contig_len_list)
    write_to_disk(nx_contig,contig_len_list, fasta_file_name, prefix)