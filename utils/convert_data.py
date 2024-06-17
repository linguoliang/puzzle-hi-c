#!/bin/env python
import pandas as pd
import sys


def convert_data(chromsize: dict, inputfilename: str, outputfilename):
    # data = pd.read_csv(chromsize, sep="\t", header=None)
    ord_scaffold = list(chromsize.keys())
    ord_scaffold.sort()
    # readname=0
    with open(inputfilename) as inputfile:
        with open(outputfilename, 'w') as outputfile:
            inputfile.readline()
            for item in inputfile:
                itemlist = item.strip().split('\t')[:8]
                if ord_scaffold.index(itemlist[1]) > ord_scaffold.index(itemlist[5]):
                    itemlist[0:4], itemlist[4:8] = itemlist[4:8], itemlist[0:4]
                # readname+=1
                outputfile.write("\t".join(itemlist) + '\n')


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("{} xx.Chrom.sizes input output!".format(sys.argv[0]))
        sys.exit(0)
    convert_data(sys.argv[1], sys.argv[2], sys.argv[3])
