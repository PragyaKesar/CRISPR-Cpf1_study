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
#import RNA
#sys.path.append("/home/titan-4/Downloads/Pragya/ViennaRNA-2.4.17/interfaces/Python3/")
data = pd.read_csv("D://PhD_related/model_evaluation/comparison_deepcpf1-crisprDT/other_datasets/HT_2/H2.csv")
data.columns = data.columns.str.strip()
#print(data.columns.tolist())
#extract columns from dataset
h1=data['Off-target']
h2=data['Query']

#define function for one-hot encoding of DNA sequence
def seq_encode():
	global dna_encode_all, dna_encodei_all, sequences
	sequences=[]
	dna_encode_all = []
	dna_encodei_all=[]
	#iterate over the columns values
	for seqs, seqs1 in zip(h1,h2):
		dna_encode = np.array([])
		dna_encodei = np.array([])
		sequences.append(seqs)
		for j in seqs:
			if j == "A":
				dna_encode = np.append(dna_encode, np.array([0,0,0,1]).T)
			if j == "T":
				dna_encode = np.append(dna_encode, np.array([0,0,1,0]).T)
			if j == "G":
				dna_encode = np.append(dna_encode, np.array([0,1,0,0]).T)
			if j == "C":
				dna_encode = np.append(dna_encode, np.array([1,0,0,0]).T)
			if j == "N":
				dna_encode = np.append(dna_encode, np.array([1,1,1,1]).T)
			if j == "-":
				dna_encode = np.append(dna_encode, np.array([0,0,0,0]).T)
		#list contains encoding of all predicted off-target sequences
		dna_encode_all.append(dna_encode)

		for j in seqs1:
			if j == "A":
				dna_encodei = np.append(dna_encodei, np.array([0,0,0,1]).T)
			if j == "T":
				dna_encodei = np.append(dna_encodei, np.array([0,0,1,0]).T)
			if j == "G":
				dna_encodei = np.append(dna_encodei, np.array([0,1,0,0]).T)
			if j == "C":
				dna_encodei = np.append(dna_encodei, np.array([1,0,0,0]).T)
			if j == "N":
				dna_encodei = np.append(dna_encodei, np.array([1,1,1,1]).T)
			if j == "-":
				dna_encodei = np.append(dna_encodei, np.array([0,0,0,0]).T)
		#list contains encoding of all gRNA sequences
		dna_encodei_all.append(dna_encodei)

#define function for calculation of positon-specific mismatches and number of mismatches
def mismatch():
	global mismatch_listF
	mismatch_listF=[]
	#iterate over the columns values
	for seqs, seqs1 in zip(h1,h2):
		mismatch_list=[]
		N=0
		#iterate over the nucleotides of the sequences
		for char1, char2 in zip(seqs, seqs1):
			if char1 == char2:
				#append 0 to the list if the nucleotide of off-target and gRNA are same
				mismatch_list=np.append(mismatch_list,np.array([0])).T
			elif char1 != char2:
				N=N+1
				#append 1 to the list if the nucleotide of off-target and gRNA are same
				mismatch_list=np.append(mismatch_list,np.array([1])).T
		mismatch_list=np.append(mismatch_list, N).T
		#list contains position-specifc mismatches of all sequence pairs
		mismatch_listF.append(mismatch_list)

#define function for calculation of presence of positon-specific bulges in gRNA and off-target sequences
def bulges():
	global bulges_OTlist, bulges_gRNAlist
	bulges_OTlist=[]
	bulges_gRNAlist=[]
	#iterate over the columns values
	for seqs, seqs1 in zip(h1,h2):
	    bulges_OT=[]
	    bulges_gRNA = []
	    n=0
	    #iterate over the nucleotides of the sequences
	    for char1, char2 in zip(seqs, seqs1):
	        if char1 == "-":
	            n=n+1
	        if char1 == "-":
	            bulges_OT=np.append(bulges_OT,np.array([1])).T
	        elif char1 != "-":
	            bulges_OT=np.append(bulges_OT,np.array([0])).T
	        if char2 == "-":
	            bulges_gRNA=np.append(bulges_gRNA,np.array([1])).T
	        elif char2 != "-":
	            bulges_gRNA=np.append(bulges_gRNA,np.array([0])).T
	    bulges_OT=np.append(bulges_OT, n).T
	    #list contains position-specifc bulges of all gRNA and off-target sequences
	    bulges_OTlist.append(bulges_OT)
	    bulges_gRNAlist.append(bulges_gRNA)

#define function for evaluating the presence of positon-specific mismatch type in gRNA and off-target sequences
def mismatch_type():
	global mm_type_list
	mm_type_list=[]
	#iterate over the columns values
	for seqs, seqs1 in zip(h1,h2):
		mm_type=[]
		for char1, char2 in zip(seqs, seqs1):
			if char1 == char2:
				mm_type=np.append(mm_type, np.array([0,0,0,0,0,0,0,0,0,0,0,0,0])).T
			elif (char1 != char2 and char1 == "A" and char2 == "T"):
				mm_type=np.append(mm_type,np.array([1,0,0,0,0,0,0,0,0,0,0,0,0])).T
			elif (char1 != char2 and char1 == "A" and char2 == "C"):
				mm_type=np.append(mm_type, np.array([0,1,0,0,0,0,0,0,0,0,0,0,0])).T
			elif(char1 != char2 and char1 == "A" and char2 == "G"):
				mm_type=np.append(mm_type,np.array([0,0,1,0,0,0,0,0,0,0,0,0,0])).T
			elif (char1 != char2 and char1 == "T" and char2 == "C"):
				mm_type=np.append(mm_type,np.array([0,0,0,1,0,0,0,0,0,0,0,0,0])).T
			elif (char1 != char2 and char1 == "T" and char2 == "G"):
				mm_type=np.append(mm_type,np.array([0,0,0,0,1,0,0,0,0,0,0,0,0])).T
			elif (char1 != char2 and char1 == "T" and char2 == "A"):
				mm_type=np.append(mm_type,np.array([0,0,0,0,0,1,0,0,0,0,0,0,0])).T
			elif (char1 != char2 and char1 == "G" and char2 == "A"):
				mm_type=np.append(mm_type,np.array([0,0,0,0,0,0,1,0,0,0,0,0,0])).T
			elif (char1 != char2 and char1 == "G" and char2 == "T"):
				mm_type=np.append(mm_type,np.array([0,0,0,0,0,0,0,1,0,0,0,0,0])).T
			elif (char1 != char2 and char1 == "G" and char2 == "C"):
				mm_type=np.append(mm_type,np.array([0,0,0,0,0,0,0,0,1,0,0,0,0])).T
			elif (char1 != char2 and char1 == "C" and char2 == "A"):
				mm_type=np.append(mm_type,np.array([0,0,0,0,0,0,0,0,0,1,0,0,0])).T
			elif (char1 != char2 and char1 == "C" and char2 == "T"):
				mm_type=np.append(mm_type,np.array([0,0,0,0,0,0,0,0,0,0,1,0,0])).T
			elif (char1 != char2 and char1 == "C" and char2 == "G"):
				mm_type=np.append(mm_type,np.array([0,0,0,0,0,0,0,0,0,0,0,1,0])).T
			elif (char1 != char2 and char1 == "N" or "-"):
				mm_type=np.append(mm_type,np.array([0,0,0,0,0,0,0,0,0,0,0,0,1])).T
			elif (char1 != char2 and char2 == "N" or "-"):
				mm_type=np.append(mm_type,np.array([0,0,0,0,0,0,0,0,0,0,0,0,1])).T
		#list contains position-specifc mismatch types between gRNA and off-target sequences
		mm_type_list.append(mm_type)

#define function for evaluating GC content of off-target sequences
def GC_contentFn():
	global GC_content_list
	GC_content_list=[]
	#iterate over the columns values
	for seqs, seqs1 in zip(h1,h2):
		length = len(seqs)
		#counts number of C in off-targt sequences
		c = seqs.count("C")
		#counts number of C in off-targt sequences
		g = seqs.count("G")
		GC = g+c
		GC_content = GC/length
		#list contains GC content of off-target sequences
		GC_content_list.append(GC_content)

#define function for evaluating dinucleotide frequencies of off-target sequences
def DiNuclFreq():
	global dinucl_freq_list1
	dinucl_freq_list1=[]
	#iterate over the columns values
	for seqs, seqs1 in zip(h1,h2):
		dinucl_freq_list=[]
		length = len(seqs)
		#counts number of different dinucleotides in off-target sequences
		AA = seqs.count("AA")
		AT = seqs.count("AT")
		AG = seqs.count("AG")
		AC = seqs.count("AC")
		TA = seqs.count("TA")
		TT = seqs.count("TT")
		TG = seqs.count("TG")
		TC = seqs.count("TC")
		CA = seqs.count("CA")
		CT = seqs.count("CT")
		CG = seqs.count("CG")
		CC = seqs.count("CC")
		GA = seqs.count("GA")
		GT = seqs.count("GT")
		GG = seqs.count("GG")
		GC = seqs.count("GC")
		dinucl_list=[AA, AT, AG, AC, TA, TT, TG, TC, CA, CT, CG, CC, GA, GT, GG, GC]
		#iterate over different dinucleotides
		for dinucl in dinucl_list:
			dinucl_freq = (dinucl/length)*100
			dinucl_freq_list.append(dinucl_freq)
		#list contains all the dinucleotide frequences of off-target sequences
		dinucl_freq_list1.append(dinucl_freq_list)

#define function for evaluating mononucleotides counts of off-target sequences
def MononuclCount():
	global mononucl_freq_list
	mononucl_freq_list=[]
	#iterate over the column values
	for seqs, seqs1 in zip(h1,h2):
		mononucl_freq_ATGC=[]
		#length of the off-target sequence
		length = len(seqs)
		#counts number of A T, G, C
		A = seqs.count("A")
		T = seqs.count("T")
		G = seqs.count("G")
		C = seqs.count("C")
		mononucl_list=[A, T, G, C]
		#iterate over different mononucleotides
		for mononucl in mononucl_list:
			mononucl_freq = (mononucl/length)*100
			mononucl_freq_ATGC.append(mononucl_freq)
		#list contains all the mononucleotide counts of off-target sequences
		mononucl_freq_list.append(mononucl_freq_ATGC)

#define function for evaluating the presence of AAAA, TTTT, GGGG, CCCC in off-target sequences
def RepetitiveSeq():
	global TetraNucl_stored_list
	TetraNucl_stored_list=[]
	#iterate over the column values
	for seqs, seqs1 in zip(h1,h2):
		TetraNucl_stored=[]
		TetraNucl_list=['AAAA', 'TTTT', 'GGGG', 'CCCC']
		#iterate over tetranucleotides
		for TetraNucl in TetraNucl_list:
			#append 0 to the list if iterated tetranucleotide is absent in off-target sequence
			if (seqs.find(TetraNucl) == -1):
				TetraNucl_stored=np.append(TetraNucl_stored,np.array([0])).T
			#append 1 to the list if iterated tetranucleotide is present in off-target sequence
			else:
				TetraNucl_stored=np.append(TetraNucl_stored,np.array([1])).T
		#list contains all the tetranucleotides of off-target sequences
		TetraNucl_stored_list.append(TetraNucl_stored)

#define function for evaluating the minimum free energy of off-target sequences
def MinFreeEnergy():
	global mfe_list
	mfe_list=[]
	#iterate over columns values
	for seqs, seqs1 in zip(h1,h2):  
	    fc = RNA.fold_compound(seqs)
	    (ss, mfe) = fc.mfe()
	    mfe = list(mfe)
	    #print("%6.2f" %(mfe), file=f)
	    #list contains minimum free energy of off-target sequences
		mfe_list.append(mfe)

#define function for evaluating the melting temperature of off-target sequences
def MeltingTemp():
	global Tm_gRNA, Tm_seed, Tm_NonSeed
	Tm_gRNA=[]
	Tm_seed=[]
	Tm_NonSeed=[]
	#iterate over the column values
	for seqs in h1:
	    seed= seqs[4:10]
	    non_seed= seqs[12:26]
	    #melting temperature of gRNA
	    Tm_gRNA1= mt.Tm_NN(seqs, nn_table=mt.R_DNA_NN1)
	    #melting teperature of seed region
	    Tm_seed1= mt.Tm_NN(seed, nn_table=mt.R_DNA_NN1)
	    #melting temperature of non-seed region
	    Tm_NonSeed1= mt.Tm_NN(non_seed, nn_table=mt.R_DNA_NN1)
	    #append the values of melting temperature of gRNA
	    Tm_gRNA.append(Tm_gRNA1)
	    #append the values of melting temperature of seed region
	    Tm_seed.append(Tm_seed1)
	    #append the values of melting temperature of non-seed region
	    Tm_NonSeed.append(Tm_NonSeed1)

#define function for the prediction of probability of off-targets
def ModelPred():
	Efficiency_score=[]
	Features = FinalDF.iloc[:,0:-1]
	Targets = FinalDF.iloc[:,-1]
	#load prediction model for AsCpf1 (AdaboostClassifier)
	Prediction = joblib.load('D://PhD_related/py_scripts/Final_scripts/AdaboostClassifier_AsCpf1.pkl')
	#prediction of probability for positive off-targets
	Efficiency_score=Prediction.predict_proba(Features)[:,1]
	h1_df=pd.DataFrame(sequences, columns=['Off-targets'])
	E_score_df=pd.DataFrame(Efficiency_score, columns=['Efficiency score'])
	#results in form of dataframe
	OutputDF = pd.concat([h1_df, E_score_df], axis=1)
	#print(OutputDF)

#function calling and dataframe construction
seq_encode()
dna_encode_df = pd.DataFrame(dna_encode_all, columns=['C_OTSeqPosition1', 'G_OTSeqPosition1', 'T_OTSeqPosition1', 'A_OTSeqPosition1', 'C_OTSeqPosition2', 'G_OTSeqPosition2', 'T_OTSeqPosition2', 'A_OTSeqPosition2', 'C_OTSeqPosition3', 'G_OTSeqPosition3', 'T_OTSeqPosition3', 'A_OTSeqPosition3', 'C_OTSeqPosition4', 'G_OTSeqPosition4', 'T_OTSeqPosition4', 'A_OTSeqPosition4', 'C_OTSeqPosition5', 'G_OTSeqPosition5', 'T_OTSeqPosition5', 'A_OTSeqPosition5', 'C_OTSeqPosition6', 'G_OTSeqPosition6', 'T_OTSeqPosition6', 'A_OTSeqPosition6', 'C_OTSeqPosition7', 'G_OTSeqPosition7', 'T_OTSeqPosition7', 'A_OTSeqPosition7', 'C_OTSeqPosition8', 'G_OTSeqPosition8', 'T_OTSeqPosition8', 'A_OTSeqPosition8', 'C_OTSeqPosition9', 'G_OTSeqPosition9', 'T_OTSeqPosition9', 'A_OTSeqPosition9', 'C_OTSeqPosition10', 'G_OTSeqPosition10', 'T_OTSeqPosition10', 'A_OTSeqPosition10', 'C_OTSeqPosition11', 'G_OTSeqPosition11', 'T_OTSeqPosition11', 'A_OTSeqPosition11', 'C_OTSeqPosition12', 'G_OTSeqPosition12', 'T_OTSeqPosition12', 'A_OTSeqPosition12', 'C_OTSeqPosition13', 'G_OTSeqPosition13', 'T_OTSeqPosition13', 'A_OTSeqPosition13', 'C_OTSeqPosition14', 'G_OTSeqPosition14', 'T_OTSeqPosition14', 'A_OTSeqPosition14', 'C_OTSeqPosition15', 'G_OTSeqPosition15', 'T_OTSeqPosition15', 'A_OTSeqPosition15', 'C_OTSeqPosition16', 'G_OTSeqPosition16', 'T_OTSeqPosition16', 'A_OTSeqPosition16', 'C_OTSeqPosition17', 'G_OTSeqPosition17', 'T_OTSeqPosition17', 'A_OTSeqPosition17', 'C_OTSeqPosition18', 'G_OTSeqPosition18', 'T_OTSeqPosition18', 'A_OTSeqPosition18', 'C_OTSeqPosition19', 'G_OTSeqPosition19', 'T_OTSeqPosition19', 'A_OTSeqPosition19', 'C_OTSeqPosition20', 'G_OTSeqPosition20', 'T_OTSeqPosition20', 'A_OTSeqPosition20', 'C_OTSeqPosition21', 'G_OTSeqPosition21', 'T_OTSeqPosition21', 'A_OTSeqPosition21', 'C_OTSeqPosition22', 'G_OTSeqPosition22', 'T_OTSeqPosition22', 'A_OTSeqPosition22', 'C_OTSeqPosition23', 'G_OTSeqPosition23', 'T_OTSeqPosition23', 'A_OTSeqPosition23', 'C_OTSeqPosition24', 'G_OTSeqPosition24', 'T_OTSeqPosition24', 'A_OTSeqPosition24', 'C_OTSeqPosition25', 'G_OTSeqPosition25', 'T_OTSeqPosition25', 'A_OTSeqPosition25', 'C_OTSeqPosition26', 'G_OTSeqPosition26', 'T_OTSeqPosition26', 'A_OTSeqPosition26', 'C_OTSeqPosition27', 'G_OTSeqPosition27', 'T_OTSeqPosition27', 'A_OTSeqPosition27'])
dna_encodei_df = pd.DataFrame(dna_encodei_all, columns=['C_QueryPosition1', 'G_QueryPosition1', 'T_QueryPosition1', 'A_QueryPosition1', 'C_QueryPosition2', 'G_QueryPosition2', 'T_QueryPosition2', 'A_QueryPosition2', 'C_QueryPosition3', 'G_QueryPosition3', 'T_QueryPosition3', 'A_QueryPosition3', 'C_QueryPosition4', 'G_QueryPosition4', 'T_QueryPosition4', 'A_QueryPosition4', 'C_QueryPosition5', 'G_QueryPosition5', 'T_QueryPosition5', 'A_QueryPosition5', 'C_QueryPosition6', 'G_QueryPosition6', 'T_QueryPosition6', 'A_QueryPosition6', 'C_QueryPosition7', 'G_QueryPosition7', 'T_QueryPosition7', 'A_QueryPosition7', 'C_QueryPosition8', 'G_QueryPosition8', 'T_QueryPosition8', 'A_QueryPosition8', 'C_QueryPosition9', 'G_QueryPosition9', 'T_QueryPosition9', 'A_QueryPosition9', 'C_QueryPosition10', 'G_QueryPosition10', 'T_QueryPosition10', 'A_QueryPosition10', 'C_QueryPosition11', 'G_QueryPosition11', 'T_QueryPosition11', 'A_QueryPosition11', 'C_QueryPosition12', 'G_QueryPosition12', 'T_QueryPosition12', 'A_QueryPosition12', 'C_QueryPosition13', 'G_QueryPosition13', 'T_QueryPosition13', 'A_QueryPosition13', 'C_QueryPosition14', 'G_QueryPosition14', 'T_QueryPosition14', 'A_QueryPosition14', 'C_QueryPosition15', 'G_QueryPosition15', 'T_QueryPosition15', 'A_QueryPosition15', 'C_QueryPosition16', 'G_QueryPosition16', 'T_QueryPosition16', 'A_QueryPosition16', 'C_QueryPosition17', 'G_QueryPosition17', 'T_QueryPosition17', 'A_QueryPosition17', 'C_QueryPosition18', 'G_QueryPosition18', 'T_QueryPosition18', 'A_QueryPosition18', 'C_QueryPosition19', 'G_QueryPosition19', 'T_QueryPosition19', 'A_QueryPosition19', 'C_QueryPosition20', 'G_QueryPosition20', 'T_QueryPosition20', 'A_QueryPosition20', 'C_QueryPosition21', 'G_QueryPosition21', 'T_QueryPosition21', 'A_QueryPosition21', 'C_QueryPosition22', 'G_QueryPosition22', 'T_QueryPosition22', 'A_QueryPosition22', 'C_QueryPosition23', 'G_QueryPosition23', 'T_QueryPosition23', 'A_QueryPosition23', 'C_QueryPosition24', 'G_QueryPosition24', 'T_QueryPosition24', 'A_QueryPosition24', 'C_QueryPosition25', 'G_QueryPosition25', 'T_QueryPosition25', 'A_QueryPosition25', 'C_QueryPosition26', 'G_QueryPosition26', 'T_QueryPosition26', 'A_QueryPosition26', 'C_QueryPosition27', 'G_QueryPosition27', 'T_QueryPosition27', 'A_QueryPosition27'])
mismatch()
mismatch_df = pd.DataFrame(mismatch_listF, columns=['mismatch_POS1', 'mismatch_POS2', 'mismatch_POS3', 'mismatch_POS4', 'mismatch_POS5', 'mismatch_POS6', 'mismatch_POS7', 'mismatch_POS8', 'mismatch_POS9', 'mismatch_POS10', 'mismatch_POS11', 'mismatch_POS12', 'mismatch_POS13', 'mismatch_POS14', 'mismatch_POS15', 'mismatch_POS16', 'mismatch_POS17', 'mismatch_POS18', 'mismatch_POS19', 'mismatch_POS20', 'mismatch_POS21', 'mismatch_POS22', 'mismatch_POS23', 'mismatch_POS24', 'mismatch_POS25', 'mismatch_POS26', 'mismatch_POS27', 'no_mismatches'])
bulges()
bulges_df = pd.DataFrame(bulges_OTlist, columns=['bulge_POS1', 'bulge_POS2', 'bulge_POS3', 'bulge_POS4', 'bulge_POS5', 'bulge_POS6', 'bulge_POS7', 'bulge_POS8', 'bulge_POS9', 'bulge_POS10', 'bulge_POS11', 'bulge_POS12', 'bulge_POS13', 'bulge_POS14', 'bulge_POS15', 'bulge_POS16', 'bulge_POS17', 'bulge_POS18', 'bulge_POS19', 'bulge_POS20', 'bulge_POS21', 'bulge_POS22', 'bulge_POS23', 'bulge_POS24', 'bulge_POS25', 'bulge_POS26', 'bulge_POS27', 'no_bulges'])
bulges_df1 = pd.DataFrame(bulges_gRNAlist, columns=['bulge_Query_POS1', 'bulge_Query_POS2', 'bulge_Query_POS3', 'bulge_Query_POS4', 'bulge_Query_POS5', 'bulge_Query_POS6', 'bulge_Query_POS7', 'bulge_Query_POS8', 'bulge_Query_POS9', 'bulge_Query_POS10', 'bulge_Query_POS11', 'bulge_Query_POS12', 'bulge_Query_POS13', 'bulge_Query_POS14', 'bulge_Query_POS15', 'bulge_Query_POS16', 'bulge_Query_POS17', 'bulge_Query_POS18', 'bulge_Query_POS19', 'bulge_Query_POS20', 'bulge_Query_POS21', 'bulge_Query_POS22', 'bulge_Query_POS23', 'bulge_Query_POS24', 'bulge_Query_POS25', 'bulge_Query_POS26', 'bulge_Query_POS27'])
mismatch_type()
mm_type_df = pd.DataFrame(mm_type_list, columns=['MM_type_A–T_POS1', 'MM_type_A–C_POS1', 'MM_type_A–G_POS1', 'MM_type_T–C_POS1', 'MM_type_T–G_POS1', 'MM_type_T–A_POS1', 'MM_type_G–A_POS1', 'MM_type_G–T_POS1', 'MM_type_G–C_POS1', 'MM_type_C–A_POS1', 'MM_type_C–T_POS1', 'MM_type_C–G_POS1', 'MM_type_other_POS1', 'MM_type_A–T_POS2', 'MM_type_A–C_POS2', 'MM_type_A–G_POS2', 'MM_type_T–C_POS2', 'MM_type_T–G_POS2', 'MM_type_T–A_POS2', 'MM_type_G–A_POS2', 'MM_type_G–T_POS2', 'MM_type_G–C_POS2', 'MM_type_C–A_POS2', 'MM_type_C–T_POS2', 'MM_type_C–G_POS2', 'MM_type_other_POS2', 'MM_type_A–T_POS3', 'MM_type_A–C_POS3', 'MM_type_A–G_POS3', 'MM_type_T–C_POS3', 'MM_type_T–G_POS3', 'MM_type_T–A_POS3', 'MM_type_G–A_POS3', 'MM_type_G–T_POS3', 'MM_type_G–C_POS3', 'MM_type_C–A_POS3', 'MM_type_C–T_POS3', 'MM_type_C–G_POS3', 'MM_type_other_POS3', 'MM_type_A–T_POS4', 'MM_type_A–C_POS4', 'MM_type_A–G_POS4', 'MM_type_T–C_POS4', 'MM_type_T–G_POS4', 'MM_type_T–A_POS4', 'MM_type_G–A_POS4', 'MM_type_G–T_POS4', 'MM_type_G–C_POS4', 'MM_type_C–A_POS4', 'MM_type_C–T_POS4', 'MM_type_C–G_POS4', 'MM_type_other_POS4', 'MM_type_A–T_POS5', 'MM_type_A–C_POS5', 'MM_type_A–G_POS5', 'MM_type_T–C_POS5', 'MM_type_T–G_POS5', 'MM_type_T–A_POS5', 'MM_type_G–A_POS5', 'MM_type_G–T_POS5', 'MM_type_G–C_POS5', 'MM_type_C–A_POS5', 'MM_type_C–T_POS5', 'MM_type_C–G_POS5', 'MM_type_other_POS5', 'MM_type_A–T_POS6', 'MM_type_A–C_POS6', 'MM_type_A–G_POS6', 'MM_type_T–C_POS6', 'MM_type_T–G_POS6', 'MM_type_T–A_POS6', 'MM_type_G–A_POS6', 'MM_type_G–T_POS6', 'MM_type_G–C_POS6', 'MM_type_C–A_POS6', 'MM_type_C–T_POS6', 'MM_type_C–G_POS6', 'MM_type_other_POS6', 'MM_type_A–T_POS7', 'MM_type_A–C_POS7', 'MM_type_A–G_POS7', 'MM_type_T–C_POS7', 'MM_type_T–G_POS7', 'MM_type_T–A_POS7', 'MM_type_G–A_POS7', 'MM_type_G–T_POS7', 'MM_type_G–C_POS7', 'MM_type_C–A_POS7', 'MM_type_C–T_POS7', 'MM_type_C–G_POS7', 'MM_type_other_POS7', 'MM_type_A–T_POS8', 'MM_type_A–C_POS8', 'MM_type_A–G_POS8', 'MM_type_T–C_POS8', 'MM_type_T–G_POS8', 'MM_type_T–A_POS8', 'MM_type_G–A_POS8', 'MM_type_G–T_POS8', 'MM_type_G–C_POS8', 'MM_type_C–A_POS8', 'MM_type_C–T_POS8', 'MM_type_C–G_POS8', 'MM_type_other_POS8', 'MM_type_A–T_POS9', 'MM_type_A–C_POS9', 'MM_type_A–G_POS9', 'MM_type_T–C_POS9', 'MM_type_T–G_POS9', 'MM_type_T–A_POS9', 'MM_type_G–A_POS9', 'MM_type_G–T_POS9', 'MM_type_G–C_POS9', 'MM_type_C–A_POS9', 'MM_type_C–T_POS9', 'MM_type_C–G_POS9', 'MM_type_other_POS9', 'MM_type_A–T_POS10', 'MM_type_A–C_POS10', 'MM_type_A–G_POS10', 'MM_type_T–C_POS10', 'MM_type_T–G_POS10', 'MM_type_T–A_POS10', 'MM_type_G–A_POS10', 'MM_type_G–T_POS10', 'MM_type_G–C_POS10', 'MM_type_C–A_POS10', 'MM_type_C–T_POS10', 'MM_type_C–G_POS10', 'MM_type_other_POS10', 'MM_type_A–T_POS11', 'MM_type_A–C_POS11', 'MM_type_A–G_POS11', 'MM_type_T–C_POS11', 'MM_type_T–G_POS11', 'MM_type_T–A_POS11', 'MM_type_G–A_POS11', 'MM_type_G–T_POS11', 'MM_type_G–C_POS11', 'MM_type_C–A_POS11', 'MM_type_C–T_POS11', 'MM_type_C–G_POS11', 'MM_type_other_POS11', 'MM_type_A–T_POS12', 'MM_type_A–C_POS12', 'MM_type_A–G_POS12', 'MM_type_T–C_POS12', 'MM_type_T–G_POS12', 'MM_type_T–A_POS12', 'MM_type_G–A_POS12', 'MM_type_G–T_POS12', 'MM_type_G–C_POS12', 'MM_type_C–A_POS12', 'MM_type_C–T_POS12', 'MM_type_C–G_POS12', 'MM_type_other_POS12', 'MM_type_A–T_POS13', 'MM_type_A–C_POS13', 'MM_type_A–G_POS13', 'MM_type_T–C_POS13', 'MM_type_T–G_POS13', 'MM_type_T–A_POS13', 'MM_type_G–A_POS13', 'MM_type_G–T_POS13', 'MM_type_G–C_POS13', 'MM_type_C–A_POS13', 'MM_type_C–T_POS13', 'MM_type_C–G_POS13', 'MM_type_other_POS13', 'MM_type_A–T_POS14', 'MM_type_A–C_POS14', 'MM_type_A–G_POS14', 'MM_type_T–C_POS14', 'MM_type_T–G_POS14', 'MM_type_T–A_POS14', 'MM_type_G–A_POS14', 'MM_type_G–T_POS14', 'MM_type_G–C_POS14', 'MM_type_C–A_POS14', 'MM_type_C–T_POS14', 'MM_type_C–G_POS14', 'MM_type_other_POS14', 'MM_type_A–T_POS15', 'MM_type_A–C_POS15', 'MM_type_A–G_POS15', 'MM_type_T–C_POS15', 'MM_type_T–G_POS15', 'MM_type_T–A_POS15', 'MM_type_G–A_POS15', 'MM_type_G–T_POS15', 'MM_type_G–C_POS15', 'MM_type_C–A_POS15', 'MM_type_C–T_POS15', 'MM_type_C–G_POS15', 'MM_type_other_POS15', 'MM_type_A–T_POS16', 'MM_type_A–C_POS16', 'MM_type_A–G_POS16', 'MM_type_T–C_POS16', 'MM_type_T–G_POS16', 'MM_type_T–A_POS16', 'MM_type_G–A_POS16', 'MM_type_G–T_POS16', 'MM_type_G–C_POS16', 'MM_type_C–A_POS16', 'MM_type_C–T_POS16', 'MM_type_C–G_POS16', 'MM_type_other_POS16', 'MM_type_A–T_POS17', 'MM_type_A–C_POS17', 'MM_type_A–G_POS17', 'MM_type_T–C_POS17', 'MM_type_T–G_POS17', 'MM_type_T–A_POS17', 'MM_type_G–A_POS17', 'MM_type_G–T_POS17', 'MM_type_G–C_POS17', 'MM_type_C–A_POS17', 'MM_type_C–T_POS17', 'MM_type_C–G_POS17', 'MM_type_other_POS17', 'MM_type_A–T_POS18', 'MM_type_A–C_POS18', 'MM_type_A–G_POS18', 'MM_type_T–C_POS18', 'MM_type_T–G_POS18', 'MM_type_T–A_POS18', 'MM_type_G–A_POS18', 'MM_type_G–T_POS18', 'MM_type_G–C_POS18', 'MM_type_C–A_POS18', 'MM_type_C–T_POS18', 'MM_type_C–G_POS18', 'MM_type_other_POS18', 'MM_type_A–T_POS19', 'MM_type_A–C_POS19', 'MM_type_A–G_POS19', 'MM_type_T–C_POS19', 'MM_type_T–G_POS19', 'MM_type_T–A_POS19', 'MM_type_G–A_POS19', 'MM_type_G–T_POS19', 'MM_type_G–C_POS19', 'MM_type_C–A_POS19', 'MM_type_C–T_POS19', 'MM_type_C–G_POS19', 'MM_type_other_POS19', 'MM_type_A–T_POS20', 'MM_type_A–C_POS20', 'MM_type_A–G_POS20', 'MM_type_T–C_POS20', 'MM_type_T–G_POS20', 'MM_type_T–A_POS20', 'MM_type_G–A_POS20', 'MM_type_G–T_POS20', 'MM_type_G–C_POS20', 'MM_type_C–A_POS20', 'MM_type_C–T_POS20', 'MM_type_C–G_POS20', 'MM_type_other_POS20', 'MM_type_A–T_POS21', 'MM_type_A–C_POS21', 'MM_type_A–G_POS21', 'MM_type_T–C_POS21', 'MM_type_T–G_POS21', 'MM_type_T–A_POS21', 'MM_type_G–A_POS21', 'MM_type_G–T_POS21', 'MM_type_G–C_POS21', 'MM_type_C–A_POS21', 'MM_type_C–T_POS21', 'MM_type_C–G_POS21', 'MM_type_other_POS21', 'MM_type_A–T_POS22', 'MM_type_A–C_POS22', 'MM_type_A–G_POS22', 'MM_type_T–C_POS22', 'MM_type_T–G_POS22', 'MM_type_T–A_POS22', 'MM_type_G–A_POS22', 'MM_type_G–T_POS22', 'MM_type_G–C_POS22', 'MM_type_C–A_POS22', 'MM_type_C–T_POS22', 'MM_type_C–G_POS22', 'MM_type_other_POS22', 'MM_type_A–T_POS23', 'MM_type_A–C_POS23', 'MM_type_A–G_POS23', 'MM_type_T–C_POS23', 'MM_type_T–G_POS23', 'MM_type_T–A_POS23', 'MM_type_G–A_POS23', 'MM_type_G–T_POS23', 'MM_type_G–C_POS23', 'MM_type_C–A_POS23', 'MM_type_C–T_POS23', 'MM_type_C–G_POS23', 'MM_type_other_POS23', 'MM_type_A–T_POS24', 'MM_type_A–C_POS24', 'MM_type_A–G_POS24', 'MM_type_T–C_POS24', 'MM_type_T–G_POS24', 'MM_type_T–A_POS24', 'MM_type_G–A_POS24', 'MM_type_G–T_POS24', 'MM_type_G–C_POS24', 'MM_type_C–A_POS24', 'MM_type_C–T_POS24', 'MM_type_C–G_POS24', 'MM_type_other_POS24', 'MM_type_A–T_POS25', 'MM_type_A–C_POS25', 'MM_type_A–G_POS25', 'MM_type_T–C_POS25', 'MM_type_T–G_POS25', 'MM_type_T–A_POS25', 'MM_type_G–A_POS25', 'MM_type_G–T_POS25', 'MM_type_G–C_POS25', 'MM_type_C–A_POS25', 'MM_type_C–T_POS25', 'MM_type_C–G_POS25', 'MM_type_other_POS25', 'MM_type_A–T_POS26', 'MM_type_A–C_POS26', 'MM_type_A–G_POS26', 'MM_type_T–C_POS26', 'MM_type_T–G_POS26', 'MM_type_T–A_POS26', 'MM_type_G–A_POS26', 'MM_type_G–T_POS26', 'MM_type_G–C_POS26', 'MM_type_C–A_POS26', 'MM_type_C–T_POS26', 'MM_type_C–G_POS26', 'MM_type_other_POS26', 'MM_type_A–T_POS27', 'MM_type_A–C_POS27', 'MM_type_A–G_POS27', 'MM_type_T–C_POS27', 'MM_type_T–G_POS27', 'MM_type_T–A_POS27', 'MM_type_G–A_POS27', 'MM_type_G–T_POS27', 'MM_type_G–C_POS27', 'MM_type_C–A_POS27', 'MM_type_C–T_POS27', 'MM_type_C–G_POS27', 'MM_type_other_POS27'])
GC_contentFn()
GC_Content_df = pd.DataFrame(GC_content_list, columns=['OT-GC_content'])
DiNuclFreq()
dinucl_df = pd.DataFrame(dinucl_freq_list1, columns=['AA_freq', 'AT_freq', 'AG_freq', 'AC_freq', 'CA_freq', 'CT_freq', 'CG_freq', 'CC_freq', 'GA_freq', 'GT_freq', 'GG_freq', 'GC_freq', 'TA_freq', 'TT_freq', 'TG_freq', 'TC_freq'])
MononuclCount()
mononucl_df = pd.DataFrame(mononucl_freq_list, columns=['No_of_A', 'No_of_T', 'No_of_G', 'No_of_C'])
RepetitiveSeq()
Tetranucl_df = pd.DataFrame(TetraNucl_stored_list, columns=['AAAA_OT', 'TTTT_OT', 'GGGG_OT', 'CCCC_OT'])
#MinFreeEnergy()
mfe_df = pd.DataFrame(mfe_list, columns=['mfe'])
MeltingTemp()
Tm_df = pd.DataFrame(Tm_gRNA, columns=['tm_gRNA'])
Tm_df1 = pd.DataFrame(Tm_seed, columns=['tm_seed'])
Tm_df2 = pd.DataFrame(Tm_NonSeed, columns=['tm_non-seed'])
Y_df=pd.DataFrame(columns=['Y'])
#merging of dataframe to construct final dataframe for prediction
FinalDF = pd.concat([dna_encode_df, dna_encodei_df, mismatch_df, bulges_df, bulges_df1, mm_type_df, GC_Content_df, dinucl_df, mononucl_df, Tetranucl_df, mfe_df, Tm_df, Tm_df1, Tm_df2, Y_df], axis=1)
pd.set_option('display.max_columns', None)

#removal of columns not present in training dataset
FinalDF.drop(['MM_type_A–C_POS1', 'MM_type_A–G_POS1', 'MM_type_T–C_POS1', 'MM_type_T–G_POS1', 'MM_type_T–A_POS1', 'MM_type_G–A_POS1', 'MM_type_G–C_POS1', 'MM_type_C–A_POS1', 'MM_type_C–G_POS1', 'MM_type_other_POS1', 'MM_type_T–A_POS27', 'MM_type_G–A_POS27', 'MM_type_G–T_POS27', 'MM_type_C–A_POS27', 'MM_type_other_POS27', 'bulge_POS1', 'bulge_POS2', 'bulge_POS3', 'bulge_POS25', 'bulge_POS26', 'bulge_POS27', 'bulge_Query_POS1', 'bulge_Query_POS2', 'bulge_Query_POS3', 'bulge_Query_POS25', 'bulge_Query_POS26', 'bulge_Query_POS27'], axis=1, inplace=True)