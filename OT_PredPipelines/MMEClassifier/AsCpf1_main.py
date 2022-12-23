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
from MMEFeatureExtractor_AsCpf1 import AsCpf1_MMEFeatureExtractor
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import sys
import RNA

#Probability prediction using AdabostClassifier with pre-defined function
AsCpf1FE = AsCpf1_MMEFeatureExtractor()
AsCpf1FE.PAMEnergy()
AsCpf1FE.MismatchEnergy()
AsCpf1FE.DataFrameConstruct()
AsCpf1FE.ModelPred()
