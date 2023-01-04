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
from MMEFeatureExtractor_LbCpf1 import LbCpf1_MMEFeatureExtractor
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import sys
import RNA

#Probability prediction using AdabostClassifier with pre-defined function
LbCpf1FE = AsCpf1_MMEFeatureExtractor()
LbCpf1FE.PAMEnergy()
LbCpf1FE.MismatchEnergy()
LbCpf1FE.DataFrameConstruct()
LbCpf1FE.ModelPred()
