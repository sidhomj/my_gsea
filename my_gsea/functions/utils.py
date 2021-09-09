
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


def GSEA(gene_sel,gene_score,p_thresh=0.05,num_perm=10,num=None):
    # # Get Enriched genes/gene sets in clusters
    # gene_set = gp.parser.gsea_gmt_parser('/home/jsidhom1/DeepSC/Data/gene_sets/h.all.v7.1.symbols.gmt')
    # gene_set = gp.parser.gsea_gmt_parser('/home/jsidhom1/DeepSC/Data/gene_sets/c5.all.v7.1.symbols.gmt')
    # gene_set =  gp.parser.gsea_gmt_parser('/home/jsidhom1/DeepSC/Data/gene_sets/c7.all.v7.1.symbols.gmt')
    # # gene_set = gp.parser.gsea_gmt_parser('/home/jsidhom1/DeepSC/Data/gene_sets/c2.cp.biocarta.v7.1.symbols.gmt')
    #
    gene_sets = ['c2.cp.biocarta.v7.4.symbols.gmt',
                 'c2.cp.kegg.v7.4.symbols.gmt',
                 'c2.cp.pid.v7.4.symbols.gmt',
                 'c2.cp.reactome.v7.4.symbols.gmt']
    #
    # gene_sets = ['c2.cp.biocarta.v7.4.symbols.gmt',
    #              'c2.cp.pid.v7.4.symbols.gmt']
    #gene_sets = ['c2.cp.kegg.v7.1.symbols.gmt']
    # gene_sets = ['c2.cp.pid.v7.1.symbols.gmt']
    # gene_sets = ['h.all.v7.4.symbols.gmt']

    gene_sets = ['c7.all.v7.4.symbols.gmt']

    gene_set = {}
    for gs in gene_sets:
        gene_set_t = gp.parser.gsea_gmt_parser(os.path.join('../Data/gene_sets',gs))
        gene_set.update(gene_set_t)

    aucs = []
    p_vals = []
    genes_in_gs = []
    gene_set_list = []
    overlap = []
    objects = []
    for gs in gene_set:
            gene_set_list.append(gs)
            obj_out = graph_object()
            auc_i,p_i,o = compute_gs_enrichment(obj_out,gene_sel,gene_score,gene_set[gs],num_perm=num_perm)
            aucs.append(auc_i)
            p_vals.append(p_i)
            overlap.append(str(len(o))+'/'+str(len(gene_set[gs])))
            genes_in_gs.append(';'.join(o))
            objects.append(obj_out)

    df_out = pd.DataFrame()
    df_out['Term'] = gene_set_list
    df_out['Overlap'] = overlap
    df_out['AUC'] = aucs
    df_out['P-value'] = p_vals
    _,df_out['Adjusted P-value'],_,_ = multipletests(p_vals,method='fdr_bh')
    df_out['Genes'] = genes_in_gs
    df_out.sort_values(by=['Adjusted P-value','AUC'],inplace=True,ascending=[True,False])
    df_out = df_out[df_out['Adjusted P-value']<p_thresh]
    return df_out, dict(zip(gene_set_list,objects))



gene_score = np.array(df_scores['1'])
gene_sel = np.array(df_rank['1'])
enr_results, do = GSEA(gene_sel,gene_score,p_thresh=0.05,num_perm=1000)
enr_results.sort_values(by='AUC',ascending=False,inplace=True)

sel = 0
matplotlib.rcdefaults()
path_sel = enr_results['Term'].iloc[sel]
fig,ax = plt.subplots(3,1,figsize=(10,5))
x = copy.copy(do[path_sel].idx)
ax[0].vlines(np.where(x!=0),ymin=0,ymax=1)
ax[0].set_xticks([])
ax[0].set_yticks([])
ax[0].set_xlim([0,len(gene_score)])
x = np.array(range(len(do[path_sel].idx)))
y = gene_score
ax[1].plot(x[y!=0],y[y!=0])
ax[1].set_xticks([])
ax[1].axhline(y=0,color='k')
ax[1].set_xlim([0,len(gene_score)])
ax[2].plot(do[path_sel].idx_sum)
ax[2].set_xlim([0,len(gene_score)])
plt.suptitle(path_sel)
plt.subplots_adjust(hspace=0.00)
