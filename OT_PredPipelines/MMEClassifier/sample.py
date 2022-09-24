class FeatureExtractor:
    def __init__(self):
        pass
    def PAMEnergy(self, x, y):
        global dna_encode_all, dna_encodei_all, sequences
        sequences=[]
        dna_encode_all = []
        dna_encodei_all=[]
	    #iterate over the columns values
        for seqs, seqs1 in zip(h1,h2):
        	dna_encode = np.array([])
        	dna_encodei = np.array([])
        	PAM = seqs1[0:4]
        	PAM1 = seqs[0:4]
        	i=0
        	for j, k in zip(PAM, PAM1):
        		i=i+1
        		if j and k == "A" and i==1:
        			dna_encode = np.append(dna_encode, np.array([1.16]).T)
        			dna_encodei = np.append(dna_encodei, np.array([1.16]).T)
        		if j and k == "A" and i==2:
        			dna_encode = np.append(dna_encode, np.array(['inf']).T)
        			dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
        	dna_encode_all.append(dna_encode)
        	dna_encodei_all.append(dna_encodei)
        	print(dna_encodei_all)


import warnings  
import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
import sys
data = pd.read_csv("D://PhD_related/model_evaluation/comparison_deepcpf1-crisprDT/other_datasets/HT_2/HEK-lenti-dataset1.csv")
data.columns = data.columns.str.strip()
#print(data.columns.tolist())
#extract columns from dataset
h1=data['Off-target']
h2=data['Query']
MME=[0.82, 0.69, 0.89, 0.87, 0.65, 0.84, 2.12, 1.20, 0.77, 1.27, -0.05, 0.71, 2.30, 1.71, 0.73, 0.49, 0.01, 0.22, 0.13, 0.04]
ClassVar=FeatureExtractor()
ClassVar.PAMEnergy(h1, h2)