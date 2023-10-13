#!/bin/env python
import pandas as pd
import sys
if len(sys.argv)!=4:
    print("{} xx.Chrom.sizes input output!".format(sys.argv[0]))
    sys.exit(0)
data=pd.read_csv(sys.argv[1],sep="\t",header=None)
ord_scaffold=list(data[0])
ord_scaffold.sort()
readname=0
with open(sys.argv[2]) as inputfile:
    with open(sys.argv[3],'w') as outputfile:
        inputfile.readline()
        for item in inputfile:
            itemlist=item.strip().split('\t')
            if ord_scaffold.index(itemlist[1])>ord_scaffold.index(itemlist[5]):
                itemlist[0:4],itemlist[4:8]=itemlist[4:8],itemlist[0:4]
            readname+=1
            outputfile.write("\t".join(itemlist)+'\n')