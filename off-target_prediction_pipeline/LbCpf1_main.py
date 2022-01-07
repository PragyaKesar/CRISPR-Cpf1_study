import warnings  
import pandas as pd
import numpy as np
np.random.seed(198)
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
import csv
import sys
from difflib import SequenceMatcher
import LbCpf1FeatureExtractor
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import sys
import RNA


#Features value calculation using pre-defined functions
LbCpf1FeatureExtractor.mismatch_type()
LbCpf1FeatureExtractor.mismatch()
LbCpf1FeatureExtractor.GC_contentFn()
LbCpf1FeatureExtractor.DiNuclFreq()
LbCpf1FeatureExtractor.RepetitiveSeq()
LbCpf1FeatureExtractor.MononuclCount()
LbCpf1FeatureExtractor.MeltingTemp()
LbCpf1FeatureExtractor.MinFreeEnergy()

#Probability prediction using AdabostClassifier with pre-defined function
LbCpf1FeatureExtractor.ModelPred()
