import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gseapy as gp
import os
import copy
from statsmodels.stats.multitest import multipletests
from my_gsea.functions import utils
from multiprocessing import Pool
from multiprocessing import get_context

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
            auc_i, p_i, obj_out,overlap_num,overlap_genes = utils.compute_gs_enrichment(gene_sel,
                                                        gene_score,
                                                        self.gene_set[gs],
                                                        num_perm=self.num_perm)
            aucs.append(auc_i)
            p_vals.append(p_i)
            overlap.append(overlap_num)
            genes_in_gs.append(overlap_genes)
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

    def Run_parallel(self,gene_sel,gene_score,num_workers=2,p=None):
        if p is None:
            p = get_context('fork').Pool(num_workers)

        gene_sets = [self.gene_set[gs] for gs in self.gene_set]
        gene_set_list = [gs for gs in self.gene_set]
        num_ins = len(gene_sets)

        args = list(zip([gene_sel]*num_ins,
                        [gene_score]*num_ins,
                        gene_sets,
                        [self.num_perm]*num_ins
            )
        )

        out = p.starmap(utils.compute_gs_enrichment,args)

        if p is None:
            p.close()
            p.join()

        out = list(zip(*out))
        aucs = out[0]
        p_vals = out[1]
        objects = out[2]
        overlap = out[3]
        genes_in_gs = out[4]

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

