from collections import OrderedDict
import pandas as pd
import argparse
# import sys
import PuzzleHiC2JBAT as JBAT

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--agp', required=True, type=str, help='Agp file.')
# parser.add_argument('-m', '--matrix', required=True, type=str, help='The matrix file path.eg: merge_nodup.txt')
parser.add_argument("-j", "--juicer_tools",required=True,type=str,help="juicer_tools path.")
# parser.add_argument('-f', '--fasta', required=True, type=str, help='Scaffold fasta file.')
parser.add_argument("-p", '--prefix', default="sample", type=str, help='Output prefix! Default: sample.')
# parser.add_argument('-s', '--binsize', default=10000, type=int, help='The bin size. Default: 10000.')
# parser.add_argument('-t', '--cutoff', default=0.3, type=float, help='Score cutoff, 0.25-0.5 recommended. default: 0.3.')
# parser.add_argument('-i', '--init_trianglesize', default=3, type=int, help='Initial triangle size. Default: 3.')
parser.add_argument('-n', '--ncpus', default=1, type=int, help='Number of threads. Default: 1.')
# parser.add_argument("-e", "--error_correction",action="store_true",help="For error correction! Default: False.")
# parser.add_argument("-g", "--gap",default=100, type=int,help="The size of gap between scaffolds. Default: 100.")
AGP_HEADER=["Chromosome", "Start", "End", "Order", "Tag", "Contig_ID", "Contig_start",
                                         "Contig_end", "Orientation"]
def agp2assembly(agp,assembly):
    with open(assembly,"w") as f:
        Chromosome_chain = OrderedDict()
        agp_df = pd.read_csv(agp, names=AGP_HEADER,sep='\t', index_col=False)
        idx=0
        agp_contigs = agp_df[agp_df.Tag == "W"]
        for i in range(len(agp_contigs)):
            Chromosome = agp_contigs.iloc[i]["Chromosome"]
            if Chromosome not in Chromosome_chain:
                Chromosome_chain[Chromosome]=[]
            # if agp_df.iloc[i].Tag == "U":
            #     Chromosome_chain[Chromosome].append(str(contig_num+1))
            idx+=1
            Contig_ID=agp_df.iloc[i]["Contig_ID"]
            Contig_end=agp_df.iloc[i]["Contig_end"]
            Orientation=agp_df.iloc[i]["Orientation"]
            print(f'>{Contig_ID} {idx} {Contig_end}',file=f)
            if Orientation == "+":
                Orientation=""
            Chromosome_chain[Chromosome].append(Orientation+str(idx))
        for Chromosome in Chromosome_chain:
            print(" ".join(map(str,Chromosome_chain[Chromosome])),file=f)
        agp_contigs=JBAT.convert_to_supter_scaffold(agp_contigs)

if __name__ == '__main__':
    args = parser.parse_args()
    fasta_file_name = args.fasta
    prefix = args.prefix
    prefix.replace(".assembly")
    # nx_list = get_nx_list(args.nx)
    # contig_len_list = get_contig_length(fasta_file_name, True)
    # nx_contig = caculate_nx(nx_list, contig_len_list)
    # print_to_screen(nx_contig, contig_len_list)
    # write_to_disk(nx_contig, contig_len_list, fasta_file_name, prefix)
    # if len(sys.argv) != 3:
    #     print("{} input.agp output.assembly".format(sys.argv[0]))
    #     sys.exit(0)
    # agp2assembly(sys.argv[1], sys.argv[2])
    agp2assembly(args.agp,"Arabidopsis.assembly")