import subprocess
from collections import OrderedDict
import pandas as pd
import argparse
# from multiprocessing import Pool, cpu_count
# import sys
import PuzzleHiC2JBAT as JBAT
import convert_data as converscript

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--agp', required=True, type=str, help='Agp file.')
parser.add_argument('-m', '--matrix', required=True, type=str, help='The matrix file path.eg: merge_nodup.txt')
parser.add_argument("-j", "--juicer_tools",required=True,type=str,help="juicer_tools path.")
parser.add_argument("-p", '--prefix', default="sample", type=str, help='Output prefix! Default: sample.')

AGP_HEADER=["Chromosome", "Start", "End", "Order", "Tag", "Contig_ID", "Contig_start",
                                         "Contig_end", "Orientation"]
def agp2assembly(agp,assembly):
    with open(assembly,"w") as f:
        Chromosome_chain = OrderedDict()
        agp_df = pd.read_csv(agp, names=AGP_HEADER,sep='\t', index_col=False)
        idx=0
        agp_contigs = agp_df[agp_df.Tag == "W"]
        agp_contigs.iloc[:,6    ] = agp_contigs.Contig_start.astype(int)
        agp_contigs.iloc[:,7] = agp_contigs.Contig_end.astype(int)
        for i in range(len(agp_contigs)):
            Chromosome = agp_contigs.iloc[i]["Chromosome"]
            if Chromosome not in Chromosome_chain:
                Chromosome_chain[Chromosome]=[]
            # if agp_df.iloc[i].Tag == "U":
            #     Chromosome_chain[Chromosome].append(str(contig_num+1))
            idx+=1
            Contig_ID=agp_contigs.iloc[i]["Contig_ID"]
            Contig_end=agp_contigs.iloc[i]["Contig_end"]
            Orientation=agp_contigs.iloc[i]["Orientation"]
            print(f'>{Contig_ID} {idx} {Contig_end}',file=f)
            if Orientation == "+":
                Orientation=""
            Chromosome_chain[Chromosome].append(Orientation+str(idx))
        for Chromosome in Chromosome_chain:
            print(" ".join(map(str,Chromosome_chain[Chromosome])),file=f)
        return agp_contigs

if __name__ == '__main__':
    args = parser.parse_args()
    # fasta_file_name = args.fasta
    prefix = args.prefix
    prefix.replace(".assembly","")
    # nx_list = get_nx_list(args.nx)
    # contig_len_list = get_contig_length(fasta_file_name, True)
    # nx_contig = caculate_nx(nx_list, contig_len_list)
    # print_to_screen(nx_contig, contig_len_list)
    # write_to_disk(nx_contig, contig_len_list, fasta_file_name, prefix)
    # if len(sys.argv) != 3:
    #     print("{} input.agp output.assembly".format(sys.argv[0]))
    #     sys.exit(0)
    # agp2assembly(sys.argv[1], sys.argv[2])
    print(args.agp)
    agp_contigs=agp2assembly(args.agp,f"{prefix}.assembly")
    super_scaffold_agp=JBAT.convert_to_super_scaffold_agp(agp_contigs)
    fake_chrom_dict, Scaffold_dict_list, scaffold_index_dict,faker_scaffold_len_dict=JBAT.get_convert_info(super_scaffold_agp)
    JBAT.convert_contact_txt(args.matrix,Scaffold_dict_list,scaffold_index_dict,fake_chrom_dict)
    chrom_size_dict = JBAT.get_chrom_size_from_agp(super_scaffold_agp)
    with open(f"{prefix}.Chrom.sizes", 'w') as outfiles:
        chrom_keys=list(chrom_size_dict.keys())
        # chrom_keys.sort()
        for scaffold in chrom_keys:
            outfiles.write("{}\t{}\n".format(scaffold, chrom_size_dict[scaffold]))
    converscript.convert_data(chrom_size_dict, f"{args.matrix}.re", f"{prefix}_super.re")
    subprocess.run("LC_ALL=C sort -k2,2 -k6,6 {0}>{1}".format(f"{prefix}_super.re",
                                                              f"{prefix}_super.re.sort"),
                   shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    subprocess.run("{0} pre {1} {2}.hic {3}".format(args.juicer_tools,
                                                    f"{prefix}_super.re.sort", prefix,
                                                    f"{prefix}.Chrom.sizes"), shell=True, check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
