#!/bin/env python
import pandas as pd
import sys

def get_chromsize(agpfilename,outputfile):
    agp=pd.read_csv(sys.argv[1],sep='\t')
    with open(sys.argv[2],"w") as outfile:
        for i in pd.Categorical(agp.Chromosome).categories:
            temp_gap=agp[agp.Chromosome==i]
            outfile.write("{}\t{}\n".format(i,max(temp_gap.End)))

if __name__=="__main__":
    if len(sys.argv) != 3:
        print("{} agpfile outfile!".format(sys.argv[0]))
        sys.exit(0)
    get_chromsize(sys.argv[1],sys.argv[2])