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
from CombFeatureExtractor_AsCpf1 import AsCpf1_AllFeatureExtractor
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import sys
import RNA

#Probability prediction using AdabostClassifier with pre-defined function
AsCpf1FE = AsCpf1_AllFeatureExtractor()
AsCpf1FE.seq_encode()
AsCpf1FE.mismatch()
AsCpf1FE.bulges()
AsCpf1FE.mismatch_type()
AsCpf1FE.GC_contentFn()
AsCpf1FE.DiNuclFreq()
AsCpf1FE.MononuclCount()
AsCpf1FE.RepetitiveSeq()
AsCpf1FE.MinFreeEnergy()
AsCpf1FE.MeltingTemp()
AsCpf1FE.PAMEnergy()
AsCpf1FE.MismatchEnergy()
AsCpf1FE.DataFrameConstruct()
AsCpf1FE.ModelPred()
