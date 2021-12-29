import sys
import numpy as np
import pandas as pd
import os
from Bio import AlignIO
from Bio import SeqFeature
from Bio import SearchIO
from Bio.SearchIO import FastaIO
from Bio.SearchIO.FastaIO import FastaM10Parser
#Input from command line
Input=sys.argv[1]
cmd1 = '../bin/fasta36 -n -E 10000 -k 1 -r [1/-1] -f 5 -g 5 /home/titan-4/Downloads/output_files/'
cmd2=' /home/titan-4/Downloads/ncbi_dataset/ncbi_dataset/data/GCF_000001405.39/Human_genome.fna -m 10 > /home/titan-4/Downloads/output_files/fasta36_hg38-m1o.txt'
#concatenating all the variable for form a FASTA36 run command
cmd=cmd1+Input+cmd2
#Run in command line
os.system(cmd)
#open the FASTA36 result file
with open("D://PhD_related/model_evaluation/comparison_deepcpf1-crisprDT/other_datasets/HT_2/HEK-plasmid_seq_hg37-m101.txt", "r+") as script:
    #format alignment results into a csv file
    for b in FastaM10Parser(script, 'fasta-m10'):
        for alignment in b:
            for c in alignment:
                for d in c:
                    #seperate alignment for hit sequences from result file generated from FASTA36
                    file = open("D://PhD_related/model_evaluation/comparison_deepcpf1-crisprDT/other_datasets/HT_2/H1.txt", 'a')
                    file.write("%s %s %s %s %s" % (d.hit, d.hit_end, d.hit_start, d.hit_strand, c.ident_pct))
                    file.write("\n")
                    #seperate alignment for gRNA sequences from result file generated from FASTA36
                    file1 = open("D://PhD_related/model_evaluation/comparison_deepcpf1-crisprDT/other_datasets/HT_2/Q1.txt", 'a')
                    file1.write("%s" % (d.query))
                    file1.write("\n")

#expression to be searched in alignment file
search1 = 'Description: '
search2 = "Seq('"
f0 = open("D://PhD_related/model_evaluation/comparison_deepcpf1-crisprDT/other_datasets/HT_2/H1.txt", "r")
f1 = open("D://PhD_related/model_evaluation/comparison_deepcpf1-crisprDT/other_datasets/HT_2/Q1.txt", "r")
f2 = open("D://PhD_related/model_evaluation/comparison_deepcpf1-crisprDT/other_datasets/HT_2/H2.csv", "a")
length1=27

list1=[]
list2=[]
list3=[]
#iterate over lines in a generated alignment file of hit sequences
for line in f0:
    #search for a expression in an iterated line
    if search1 in line:
        line = line.replace("Description: Homo sapiens ", "")
        line = line.replace(", GRCh37.p13 Primary Assembly", "")
        chr_num= line.replace("\n", "")
        list1.append(chr_num)
    #search for another expression in an iterated line
    if search2 in line:
        line = line.replace("Seq('", "")
        line = line.replace("')", "")
        line = line.replace(" ", ", ")
        Sequence=line.replace("\n", "")
        Sequence=Sequence.split(",")
        list2.append(Sequence)
#iterate over lines in a generated alignment file of gRNA sequences
for line2 in f1:
    #search for an expression in an iterated line
    if search2 in line2:
        line2 = line2.replace("Seq('", "")
        line2 = line2.replace("')", "")
        line2=line2.replace("\n", "")
        list3.append(line2)
#display all the rows and columns
pd.set_option("display.max_rows", None, "display.max_columns", None)
#dataframe with gRNA
Query_df = pd.DataFrame(list3, columns=['gRNA'])
#dataframe with off-targets, end, start, strand, identity
Seq_df1 = pd.DataFrame(list2, columns=['Off-target', 'END', 'START', 'Strand', 'Percentage identity'])
#dataframe with chromosome number
chr_df = pd.DataFrame(list1, columns=['Chromosome number'])
#concatenate all the dataframes to form final dataframe
OT_df = pd.concat([Query_df, Seq_df1, chr_df], axis=1)
#extract specific column values
OT=OT_df['Off-target']
gRNA=OT_df['gRNA']
n=0
#iterate over the column values
for i, j in zip(OT, gRNA):
    n=n+1
    #filter sequence with length = 27
    if(len(i) == length1) and (len(i) == length1):
        pass
    else:
        #drop columns if length of gRNA and  off-target sequence is not equal to 27
        OT_df.drop(index=n)
#print final dataframe to command prompt screen
print(OT_df)