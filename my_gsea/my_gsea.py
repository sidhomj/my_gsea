import pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr
import gseapy as gp
import os
import matplotlib
import copy
from statsmodels.stats.multitest import multipletests
import scipy

class GSEA(object):
    def __init__(self,p_thresh = 0.05, num_perm = 10):
        self.p_thresh = p_thresh
        self.num_perm = num_perm

    def Load_Gene_Sets(self,gs):
        check=1

    def Run(self,gene_sel,gene_score):
        check=1


# gene_sel = np.array(DFs[cluster_sel]['genes'])
# gene_score = np.array(DFs[cluster_sel]['fc'])
# enr_results_2,do = GSEA(gene_sel,gene_score,num_perm=1000,p_thresh=0.05)
