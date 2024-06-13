from collections import OrderedDict
import pandas as pd
import sys
import utils.PuzzleHiC2JBAT as JBAT

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
        # agp_contigs=JBAT.convert_to_supter_scaffold(agp_contigs)

if __name__ == '__main__':
    # if len(sys.argv) != 3:
    #     print("{} input.agp output.assembly".format(sys.argv[0]))
    #     sys.exit(0)
    # agp2assembly(sys.argv[1], sys.argv[2])
    agp2assembly("Arabidopsis.agp","Arabidopsis.assembly")