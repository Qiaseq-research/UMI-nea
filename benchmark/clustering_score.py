import sys
from skbio.alignment import StripedSmithWaterman
from multiprocessing import Pool
from sklearn import metrics
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

def find_score(c_in0,c_in1):
    df0 = pd.read_csv(c_in0, sep=" ", header=None)
    df1 = pd.read_csv(c_in1, sep=" ", header=None)
    read_parent_id = df0[1].tolist()
    read_labels = df1[1].tolist()
    hom = metrics.homogeneity_score(read_parent_id, read_labels)
    comp = metrics.completeness_score(read_parent_id, read_labels)
    v = metrics.v_measure_score(read_parent_id, read_labels)
    print(c_in1,"homogeneity_score is", hom,"completeness_score is",comp, "V is",v)

f1 = sys.argv[1]
f2 = sys.argv[2]
find_score(f1,f2)