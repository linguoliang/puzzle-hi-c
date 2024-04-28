import sys

from Bio import SeqIO
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def GenerateChrom(agp,Contigs_dict,Chr_name="Chr02"):
#     N_seq=Seq("N"*Seq_N_len)
    Chr=SeqRecord("",id=Chr_name,name=Chr_name,description=Chr_name)
    for item in agp.values:
        # 判断Tag W 为序列 U为 N
        if item[4]=="W":
            if item[5] in Contigs_dict:
                if len(Contigs_dict[item[5]])!=item[-2]:
                    print("error!{} {}".format(Chr_name,Contigs_dict[item[5]].name))
                if item[-1]==0:
                    Chr.seq+=Contigs_dict[item[5]].seq
                else:
                    Chr.seq+=Contigs_dict[item[5]].seq.reverse_complement()
                # Chr.seq+=Seq("N"*gap)
            else:
                print("Agp contains unkown contig {}!".format(item[5]))
        else:
            Chr.seq+=Seq("N"*int(item[-4]))
    return Chr
def get_scaffold_seq(filename):
    scaffold_dict={}
    Scaffold_level_file=SeqIO.parse(filename,"fasta")
    for seq in Scaffold_level_file:
        if seq.name in scaffold_dict:
            print("error! same name")
            print(scaffold_dict[scaffold_dict[seq.name]])
            print(seq)
        else:
            scaffold_dict[seq.name]=seq
    return scaffold_dict

def main(Path,seq_data_path,res):
    all_agp = pd.read_csv(Path, sep='\t', index_col=0)
    # all_agp.to_csv("{}/Puzzle_{}.agp".format(result_path, res), sep='\t', index=False)
    scaffold_dict = get_scaffold_seq(seq_data_path)
    chroms = pd.Categorical(all_agp.Chromosome).categories
    chrom_keys = list(chroms)
    chrom_keys.sort(key=lambda x: int(x[9:]))
    chromlist = []
    for chrom in chrom_keys:
        tmpagp = all_agp[all_agp.Chromosome == chrom]
        chromlist.append(GenerateChrom(tmpagp, scaffold_dict, chrom))
    SeqIO.write(chromlist, "{}.fa".format(res), "fasta")

if __name__=="__main__":
    Path = sys.argv[1]
    seq_data_path = sys.argv[2]
    if len(sys.argv)==4:
        res=sys.argv[3]
    else:
        res="puzzle_hic"
    main(Path,seq_data_path,res)