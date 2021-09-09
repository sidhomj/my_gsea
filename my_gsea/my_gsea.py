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
        self.DO = dict(zip(gene_set_list, objects))
        self.enr_results = df_out
        self.gene_score = gene_score

    def Vis(self,path_sel):
        fig, ax = plt.subplots(3, 1, figsize=(10, 5))
        x = copy.copy(self.DO[path_sel].idx)
        ax[0].vlines(np.where(x != 0), ymin=0, ymax=1)
        ax[0].set_xticks([])
        ax[0].set_yticks([])
        ax[0].set_xlim([0, len(self.gene_score)])
        x = np.array(range(len(self.DO[path_sel].idx)))
        y = self.gene_score
        ax[1].plot(x[y != 0], y[y != 0])
        ax[1].set_xticks([])
        ax[1].axhline(y=0, color='k')
        ax[1].set_xlim([0, len(self.gene_score)])
        ax[2].plot(self.DO[path_sel].idx_sum)
        ax[2].set_xlim([0, len(self.gene_score)])
        plt.suptitle(path_sel)
        plt.subplots_adjust(hspace=0.00)
        return fig,ax

