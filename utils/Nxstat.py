#!/usr/bin/env python
# coding: utf-8

import argparse
import sys

import numpy as np
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-p", '--prefix', default="Nx", type=str, help='Prefix! default:Nx')
parser.add_argument('-f', '--fasta', required=True, type=str, help='Scaffold fasta file!')
parser.add_argument('-N', '--nx', default="50,90,95,98,99", type=str, help='Nx list! default:50,90,95,98,99')


def get_nx_list(string):
    for c in string:
        if c not in "0123456789,":
            print(f"Nx list contains invalid character: {c}!")
            sys.exit(-1)
    clist = string.split(',')
    Nxlist = [int(x) for x in clist]
    return Nxlist


def get_contig_length(fasta_file_name: str, sorted: bool = False):
    fastaIO = SeqIO.parse(fasta_file_name, "fasta")
    contig_len_list = []
    for seq in fastaIO:
        contig_len_list.append(len(seq))
    if sorted:
        contig_len_list.sort(reverse=True)
    return contig_len_list


def caculate_nx(nx_list, contig_len_list):
    sumN = sum(contig_len_list)
    len_array = np.array(contig_len_list)
    len_cumulate = len_array.cumsum().astype(float)
    percent = len_cumulate * 100 / sumN
    nx_contig = []
    contig_num = len(len_array)
    for Nx in nx_list:
        index = sum(percent < Nx)
        if index >= contig_num:
            index = contig_num - 1
        nx_contig.append([Nx, len_array[index],index+1])
    return nx_contig


def print_to_screen(nx_contig, nx_length):
    print(f"# Total size:{sum(nx_length)}\t Total scaffolds: {len(nx_length)}\n")
    for Nx, contig_len,L in nx_contig:
        print(f'N{Nx} is {contig_len}, L: {L}')


def write_to_disk(nx_contig,nx_length, fasta_file_name, prefix):
    with open(f'{fasta_file_name}_{prefix}.txt', 'w') as outfile:
        outfile.write(f"# Total size:{sum(nx_length)}\t Total scaffolds: {len(nx_length)}\n")
        outfile.write("# Nx\tLength\tL\n")
        for Nx, contig_len,L in nx_contig:
            outfile.write(f'N{Nx}\t{contig_len}\t{L}\n')


if __name__ == "__main__":
    args = parser.parse_args()
    fasta_file_name = args.fasta
    prefix = args.prefix
    nx_list = get_nx_list(args.nx)
    contig_len_list = get_contig_length(fasta_file_name, True)
    nx_contig = caculate_nx(nx_list, contig_len_list)
    print_to_screen(nx_contig,contig_len_list)
    write_to_disk(nx_contig,contig_len_list, fasta_file_name, prefix)
