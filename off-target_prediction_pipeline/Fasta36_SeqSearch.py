import sys
import numpy as np
import pandas as pd
import os
#from fastai.imports import *
from Bio import AlignIO
from Bio import SeqFeature
from Bio import SearchIO
from Bio.SearchIO import FastaIO
from Bio.SearchIO.FastaIO import FastaM10Parser
#Input from command line
Input=sys.argv[1]
Input1=sys.argv[2]
Out_Input=sys.argv[3]
Fasta36_Out = os.path.dirname(os.path.abspath(Input))
cmd1 = '../bin/fasta36 -n -E 10000 -k 1 -r [1/-1] -f 5 -g 5 '
cmd2=' '
cmd3=' -m 10 > '
cmd4='\\fasta36_hg38-m1o.txt'
#concatenating all the variable for form a FASTA36 run command
cmd=cmd1+Input+cmd2+Input1+cmd3+Fasta36_Out+cmd4
print(cmd)
#Run in command line
os.system(cmd)
cmd5='\\H1.txt'
cmd6='\\Q1.txt'
file_out=Fasta36_Out+cmd4
file_out1=Fasta36_Out+cmd5
file_out2=Fasta36_Out+cmd6
#open the FASTA36 result file
with open(file_out, "r+") as script:
    #format alignment results into a csv file
    for b in FastaM10Parser(script, 'fasta-m10'):
        for alignment in b:
            for c in alignment:
                for d in c:
                    #seperate alignment for hit sequences from result file generated from FASTA36
                    file = open(file_out1, 'a')
                    file.write("%s %s %s %s %s" % (d.hit, d.hit_end, d.hit_start, d.hit_strand, c.ident_pct))
                    file.write("\n")
                    #seperate alignment for gRNA sequences from result file generated from FASTA36
                    file1 = open(file_out2, 'a')
                    file1.write("%s" % (d.query))
                    file1.write("\n")

#expression to be searched in alignment file
search1 = 'Description: '
search2 = "Seq('"
f0 = open(file_out1, "r")
f1 = open(file_out2, "r")
f2 = open(Out_Input, "a")
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
#pd.set_option("display.max_rows", None, "display.max_columns", None)
#dataframe with gRNA
Query_df = pd.DataFrame(list3, columns=['gRNA'])
#dataframe with off-targets, end, start, strand, identity
Seq_df1 = pd.DataFrame(list2, columns=['Off-target', 'END', 'START', 'Strand', 'Percentage identity'])
h1=Seq_df1["Off-target"]
h2=Query_df["gRNA"]
LenO=[]
LenG=[]
for seqs1, seqs2 in zip(h1, h2): 
    LenG.append(len(seqs2))
    LenO.append(len(seqs1))

#dataframe with length
len_df = pd.DataFrame(LenO, columns=['Off-target Length'])
len_df2 = pd.DataFrame(LenG, columns=['gRNA Length'])
#dataframe with chromosome number
chr_df = pd.DataFrame(list1, columns=['Chromosome number'])
#concatenate all the dataframes to form final dataframe
OT_df = pd.concat([Query_df, Seq_df1, len_df, len_df2, chr_df], axis=1)
#extract specific column values
OT_df1=OT_df[(OT_df['Off-target Length']) & (OT_df['gRNA Length']) == 27]
#print final dataframe to command prompt screen
print(OT_df1)
