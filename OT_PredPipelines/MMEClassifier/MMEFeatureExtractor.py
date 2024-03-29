from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import warnings  
import pandas as pd
import numpy as np
import joblib
np.random.seed(123)
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
import sys
import Fasta36_SeqSearch
import RNA
#sys.path.append("/home/titan-4/Downloads/Pragya/ViennaRNA-2.4.17/interfaces/Python3/")
#data = pd.read_csv("D://PhD_related/model_evaluation/comparison_deepcpf1-crisprDT/other_datasets/HT_2/H2.csv")
#data.columns = data.columns.str.strip()
#print(data.columns.tolist())
h1=OT_df1['Off-target']
h2=OT_df1['gRNA']

#define function for one-hot encoding of DNA sequence
MME=[0.82, 0.69, 0.89, 0.87, 0.65, 0.84, 2.12, 1.20, 0.77, 1.27, -0.05, 0.71, 2.30, 1.71, 0.73, 0.49, 0.01, 0.22, 0.13, 0.04]

class FeatureExtractor:
    def __init__(self):
        pass

	def PAMEnergy():
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
				j=j+1
				if j and k == "A" and i==1:
					dna_encode = np.append(dna_encode, np.array([1.16]).T)
					dna_encodei = np.append(dna_encodei, np.array([1.16]).T)
				if j and k == "A" and i==2:
					dna_encode = np.append(dna_encode, np.array(['inf']).T)
					dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
				if j and k == "A" and i==3:
					dna_encode = np.append(dna_encode, np.array(['inf']).T)
					dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
				if j and k == "A" and i==4:
					dna_encode = np.append(dna_encode, np.array([-0.20]).T)
					dna_encodei = np.append(dna_encodei, np.array([-0.20]).T)
				if j and k == "T" and i==1:
					dna_encode = np.append(dna_encode, np.array([0.07]).T)
					dna_encodei = np.append(dna_encodei, np.array([0.07]).T)
				if j and k == "T" and i==2:
					dna_encode = np.append(dna_encode, np.array([0]).T)
					dna_encodei = np.append(dna_encodei, np.array([0]).T)
				if j and k == "T" and i==3:
					dna_encode = np.append(dna_encode, np.array([0]).T)
					dna_encodei = np.append(dna_encodei, np.array([0]).T)
				if j and k == "T" and i==4:
					dna_encode = np.append(dna_encode, np.array([2.21]).T)
					dna_encodei = np.append(dna_encodei, np.array([2.21]).T)
				if j and k == "G" and i==1:
					dna_encode = np.append(dna_encode, np.array([0.46]).T)
					dna_encodei = np.append(dna_encodei, np.array([0.46]).T)
				if j and k == "G" and i==2:
					dna_encode = np.append(dna_encode, np.array(['inf']).T)
					dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
				if j and k == "G" and i==3:
					dna_encode = np.append(dna_encode, np.array(['inf']).T)
					dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
				if j and k == "G" and i==4:
					dna_encode = np.append(dna_encode, np.array([0.13]).T)
					dna_encodei = np.append(dna_encodei, np.array([0.13]).T)
				if j and k == "C" and i==1:
					dna_encode = np.append(dna_encode, np.array([0.83]).T)
					dna_encodei = np.append(dna_encodei, np.array([0.83]).T)
				if j and k == "C" and i==2:
					dna_encode = np.append(dna_encode, np.array(['inf']).T)
					dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
				if j and k == "C" and i==3:
					dna_encode = np.append(dna_encode, np.array(['inf']).T)
					dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
				if j and k == "C" and i==4:
					dna_encode = np.append(dna_encode, np.array([0.43]).T)
					dna_encodei = np.append(dna_encodei, np.array([0.43]).T)
				if j and k == "N" and i==4:
					dna_encode = np.append(dna_encode, np.array([0]).T)
					dna_encodei = np.append(dna_encodei, np.array([0]).T)
				if j and k == "N" and i==3:
					dna_encode = np.append(dna_encode, np.array([0]).T)
					dna_encodei = np.append(dna_encodei, np.array([0]).T)
				if j and k == "N" and i==2:
					dna_encode = np.append(dna_encode, np.array([0]).T)
					dna_encodei = np.append(dna_encodei, np.array([0]).T)
				if j and k == "-" and i==1:
					dna_encode = np.append(dna_encode, np.array(['inf']).T)
					dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
				if j and k == "-" and i==2:
					dna_encode = np.append(dna_encode, np.array(['inf']).T)
					dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
				if j and k == "-" and i==3:
					dna_encode = np.append(dna_encode, np.array(['inf']).T)
					dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
				if j and k == "-" and i==4:
					dna_encode = np.append(dna_encode, np.array(['inf']).T)
					dna_encodei = np.append(dna_encodei, np.array(['inf']).T)
			PAMBindingEnergyQ=sum(dna_encodei)
			PAMBindingEnergyOT=sum(dna_encode)
			RelativePAM=PAMBindingEnergyOT-PAMBindingEnergyQ
			#list contains encoding of all predicted off-target sequences
			BindingEnergy.append(PAMBindingEnergyOT, PAMBindingEnergyQ, RelativePAM)
			dna_encode_all.append(dna_encode)
			dna_encodei_all.append(dna_encodei)

	

	def MismatchEnergy():
		global mismatch_listF
		mismatch_listF=[]
		MME=[0.82, 0.69, 0.89, 0.87, 0.65, 0.84, 2.12, 1.20, 0.77, 1.27, -0.05, 0.71, 2.30, 1.71, 0.73, 0.49, 0.01, 0.22, 0.13, 0.04]
		#iterate over the columns values

		for seqs, seqs1 in zip(h1,h2):
			seedO = seqs[4:10]
			trunkO= seqs[10:22]
			promiscuousO= seqs[22:28]
			seedQ= seqs1[4:10]
			trunkQ= seqs1[10:22]
			promiscuousQ= seqs1[22:28]
			mismatch_list_all=[]
			mismatch_listS=[]
			mismatch_listT=[]
			mismatch_listP=[]
			MMEvalue=0
			N=-1
			for char1, char2 in zip(seedO, seedQ):
				N=N+1
				if char1 == char2:
					#append 0 to the list if the nucleotide of off-target and gRNA are same
					mismatch_listS=np.append(mismatch_listS,np.array([0])).T
				elif char1 != char2:
					MMEvalue= MME[N]
					#print(MMEvalue)
					#append 1 to the list if the nucleotide of off-target and gRNA are same
					mismatch_listS=np.append(mismatch_listS,np.array([MMEvalue])).T
			TotalMMESeed=sum(mismatch_listS)
			N=-1
			for char1, char2 in zip(TrunkO, TrunkQ):
				N=N+1
				if char1 == char2:
					#append 0 to the list if the nucleotide of off-target and gRNA are same
					mismatch_listT=np.append(mismatch_listT,np.array([0])).T
				elif char1 != char2:
					MMEvalue= MME[N]
					#print(MMEvalue)
					#append 1 to the list if the nucleotide of off-target and gRNA are same
					mismatch_listT=np.append(mismatch_listT,np.array([MMEvalue])).T
			TotalMMETrunk=sum(mismatch_listT)
			N=-1
			for char1, char2 in zip(promiscuousO, promiscuousQ):
				N=N+1
				if char1 == char2:
					#append 0 to the list if the nucleotide of off-target and gRNA are same
					mismatch_listP=np.append(mismatch_listP,np.array([0])).T
				elif char1 != char2:
					MMEvalue= MME[N]
					#print(MMEvalue)
					#append 1 to the list if the nucleotide of off-target and gRNA are same
					mismatch_listP=np.append(mismatch_listP,np.array([MMEvalue])).T
			TotalMMEPromiscuous=sum(mismatch_listP)
		mismatch_list_all=np.append(mismatch_listS,mismatch_listT, mismatch_listP, TotalMMESeed, TotalMMETrunk, TotalMMEPromiscuous).T


#function calling and dataframe construction
PAMEnergy()
PAMOT_df = pd.DataFrame(dna_encode_all, columns=['PAM_OT_POS1', 'PAM_OT_POS2', 'PAM_OT_POS3', 'PAM_OT_POS4'])
PAMQ_df = pd.DataFrame(dna_encodei_all, columns=['PAM_Q_POS1', 'PAM_Q_POS2', 'PAM_Q_POS3', 'PAM_Q_POS4'])
BindingEnergy_df = pd.DataFrame(BindingEnergy, columns=['PAM binding energy OT', 'PAM binding energy Query', 'relative PAM binding energy']
MismatchEnergy()
MME_df = pd.DataFrame(mismatch_list_all, columns=['MM_SeedR_POS1', 'MM_SeedR_POS2', 'MM_SeedR_POS3', 'MM_SeedR_POS4', 'MM_SeedR_POS5', 'MM_SeedR_POS6', 'MM_TrunkR_POS7', 'MM_TrunkR_POS8', 'MM_TrunkR_POS9', 'MM_TrunkR_POS10', 'MM_TrunkR_POS11', 'MM_TrunkR_POS12', 'MM_TrunkR_POS13', 'MM_TrunkR_POS14', 'MM_TrunkR_POS15', 'MM_TrunkR_POS16', 'MM_TrunkR_POS17', 'MM_TrunkR_POS18', 'MM_PromiscuousR_POS19', 'MM_PromiscuousR_POS20', 'MM_PromiscuousR_POS21', 'MM_PromiscuousR_POS22', 'MM_PromiscuousR_POS23', 'Total MME ssed', 'total MME trunk', 'total MME promiscuous', 'total MME gRNA'])
Y_df=pd.DataFrame(columns=['Y'])
#merging of dataframe to construct final dataframe for prediction
FinalDF = pd.concat([PAMOT, PAMQ, BindingEnergy_df, MME_df, Y_df], axis=1)
pd.set_option('display.max_columns', None)