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
from my_gsea.functions import utils

class GSEA(object):
    def __init__(self,p_thresh = 0.05, num_perm = 10):
        self.p_thresh = p_thresh
        self.num_perm = num_perm

    def Load_Gene_Sets(self,gene_sets):
        gene_set = {}
        dir_path = os.path.dirname(os.path.realpath(__file__))
        for gs in gene_sets:
            gene_set_t = gp.parser.gsea_gmt_parser(os.path.join(dir_path,'gene_sets', gs))
            gene_set.update(gene_set_t)
        self.gene_set = gene_set

    def Run(self,gene_sel,gene_score):
        aucs = []
        p_vals = []
        genes_in_gs = []
        gene_set_list = []
        overlap = []
        objects = []
        for gs in self.gene_set:
            gene_set_list.append(gs)
            obj_out = utils.graph_object()
            auc_i, p_i, o = utils.compute_gs_enrichment(obj_out,
                                                        gene_sel,
                                                        gene_score,
                                                        self.gene_set[gs],
                                                        num_perm=self.num_perm)
            aucs.append(auc_i)
            p_vals.append(p_i)
            overlap.append(str(len(o)) + '/' + str(len(self.gene_set[gs])))
            genes_in_gs.append(';'.join(o))
            objects.append(obj_out)

        df_out = pd.DataFrame()
        df_out['Term'] = gene_set_list
        df_out['Overlap'] = overlap
        df_out['AUC'] = aucs
        df_out['P-value'] = p_vals
        _, df_out['Adjusted P-value'], _, _ = multipletests(p_vals, method='fdr_bh')
        df_out['Genes'] = genes_in_gs
        df_out.sort_values(by=['Adjusted P-value', 'AUC'], inplace=True, ascending=[True, False])
        df_out = df_out[df_out['Adjusted P-value'] < self.p_thresh]
        return df_out, dict(zip(gene_set_list, objects))
        check=1


# gene_sel = np.array(DFs[cluster_sel]['genes'])
# gene_score = np.array(DFs[cluster_sel]['fc'])
# enr_results_2,do = GSEA(gene_sel,gene_score,num_perm=1000,p_thresh=0.05)
