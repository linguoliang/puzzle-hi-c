from collections import OrderedDict
import pandas as pd
import sys

AGP_HEADER=["Chromosome", "Start", "End", "Order", "Tag", "Contig_ID", "Contig_start",
                                         "Contig_end", "Orientation"]
def agp2assembly(agp,assembly):
    Chromosome_chain=OrderedDict()
    agp_df = pd.DataFrame(agp, columns=AGP_HEADER, header=None, index=None)
    idx=0
    agp_contigs = agp_df[agp_df.Tag == "W"]
    agp_gaps = agp_df[agp_df.Tag == "U"]
    contig_num = len(agp_contigs)
    gap_num=len(agp_gaps)
    if gap_num > 0:
        gap_length = agp_gaps.iloc[0,5]
    for i in range(len(agp_df)):
        Chromosome = agp_df.iloc[i]["Chromosome"]
        if Chromosome not in Chromosome_chain:
            Chromosome_chain[Chromosome]=[]
        if agp_df.iloc[i].Tag == "U":
            Chromosome_chain[Chromosome].append(str(contig_num+1))
        else:
            Contig_ID=agp_df.iloc[i]["Contig_ID"]
            Contig_end=agp_df.iloc[i]["Contig_end"]
            Orientation=agp_df.iloc[i]["Orientation"]
            print(f'>{Contig_ID} {idx} {Contig_end}',file=assembly)
            if Orientation == "+":
                Orientation=""
            Chromosome_chain[Chromosome].append(Orientation+str(idx+1))
    if gap_num > 0:
        print(f'>hic_gap_{contig_num+1} {contig_num+1} {gap_length}',file=assembly)
    for Chromosome in Chromosome_chain:
        print(" ".join(map(str,Chromosome_chain[Chromosome])),file=assembly)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("{} input.agp output.assembly".format(sys.argv[0]))
        sys.exit(0)
    agp2assembly(sys.argv[1], sys.argv[2])
    # agp2assembly(agpfile,assembly)