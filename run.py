from my_gsea.my_gsea import GSEA
import pickle
import numpy as np

with open('top_genes.pkl','rb') as f:
    df_rank,df_scores = pickle.load(f)

gene_sel = np.array(df_rank['1'])
gene_score = np.array(df_scores['1'])
gene_sets = ['h.all.v7.4.symbols.gmt',
             'c2.cp.kegg.v7.4.symbols.gmt',
             'c2.cp.biocarta.v7.4.symbols.gmt',
             'c2.cp.wikipathways.v7.4.symbols.gmt']
gene_sets = ['c7.all.v7.4.symbols.gmt']

gs_obj = GSEA()
gs_obj.Load_Gene_Sets(gene_sets)
gs_obj.Run_parallel(gene_sel,gene_score,10)
gs_obj.Vis(gs_obj.enr_results['Term'].iloc[1])