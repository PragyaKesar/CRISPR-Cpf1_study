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
import AsCpf1FeatureExtractor
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import sys
import RNA

#Features value calculation using pre-defined functions
AsCpf1FeatureExtractor.mismatch_type()
AsCpf1FeatureExtractor.mismatch()
AsCpf1FeatureExtractor.GC_contentFn()
AsCpf1FeatureExtractor.DiNuclFreq()
AsCpf1FeatureExtractor.RepetitiveSeq()
AsCpf1FeatureExtractor.MononuclCount()
AsCpf1FeatureExtractor.MeltingTemp()
AsCpf1FeatureExtractor.MinFreeEnergy()

#Probability prediction using AdabostClassifier with pre-defined function
AsCpf1FeatureExtractor.ModelPred()
