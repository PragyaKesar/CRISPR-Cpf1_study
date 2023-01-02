import warnings  
import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
import csv
import sys
from difflib import SequenceMatcher
from CombFeatureExtractor_AsCpf1 import LbCpf1_AllFeatureExtractor
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import sys
import RNA

#Probability prediction using AdabostClassifier with pre-defined function
LbCpf1FE = LbCpf1_AllFeatureExtractor()
LbCpf1FE.seq_encode()
LbCpf1FE.mismatch()
LbCpf1FE.bulges()
LbCpf1FE.mismatch_type()
LbCpf1FE.GC_contentFn()
LbCpf1FE.DiNuclFreq()
LbCpf1FE.MononuclCount()
LbCpf1FE.RepetitiveSeq()
LbCpf1FE.MinFreeEnergy()
LbCpf1FE.MeltingTemp()
LbCpf1FE.PAMEnergy()
LbCpf1FE.MismatchEnergy()
LbCpf1FE.DataFrameConstruct()
LbCpf1FE.ModelPred()
