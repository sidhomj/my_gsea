import numpy as np
import scipy

class graph_object(object):
    def __init__(self):
        check=1

def compute_gs_enrichment(obj_out,gene_sel,gene_score,set,num_perm=500):
    overlap,idx_1,idx_2 = np.intersect1d(gene_sel,set,return_indices=True)
    overlap_out = np.flip(overlap[np.argsort(gene_score[idx_1])])
    if len(overlap)>0:
        idx = np.zeros_like(gene_sel)
        idx[idx_1] = gene_score[idx_1]
        obj_out.idx = idx
        idx = np.cumsum(idx)
        obj_out.idx_sum = idx
        auc = scipy.integrate.simps(idx,range(len(idx)))#/len(overlap_out)#/idx[-1]#/(np.max(idx)*len(idx))
        auc_test = auc
        aucs = []
        for _ in range(num_perm):
            idx_1 = np.random.choice(len(gene_sel),len(idx_1),replace=False)
            idx = np.zeros_like(gene_sel)
            idx[idx_1] = gene_score[idx_1]
            idx = np.cumsum(idx)
            auc = scipy.integrate.simps(idx, range(len(idx)))#/len(overlap_out)#/idx[-1] #/ (np.max(idx) * len(idx))
            aucs.append(auc)
        return auc_test, (1-scipy.stats.percentileofscore(aucs,auc_test)/100),overlap_out
    else:
        return 0.5, 1.0, []


