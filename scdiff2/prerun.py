#!/usr/bin/env python
# coding: utf-8
# Author: Jun Ding
# Email: junding (at) cs (dot) cmu (dot) edu
# Date: June. 29th, 2020
# 
# This scdiff software suite is desinged to infer the clusters, trajectories, and regulatory
# networks underlying dynamic biological process (e.g., cell differntiation, disease progression)
# based on given time-series single-cell expression input data.  Please use "scdiff -h" for the detailed usage.
# 
# This scdiff prerun program use scanpy package to learn the initial clusters/trajectories, which will be used as the input  
# to the scdiff2 main program to learn the detailed underlying regulatory networks and refined trajectories. 
#
# This software is freely avaible for academic uses. 
# For any commerical usage, please contact me at the email address above.
# All rights reserved.
# Please don NOT modify the above statement.




# In[1]:
import pdb,sys,os
import anndata
import scanpy as sc
from File import *
import pandas as pd
import argparse
import matplotlib
matplotlib.use('Agg')

def prerun(exFn,outdir,iformat,mindisp,cluRes):
    # # read in tab.txt file and save it to h5file 
    if os.path.exists(outdir)==False:
        os.mkdir(outdir)

    TabFile(exFn).toH5("\t","%s/%s"%(outdir,exFn.split("/")[-1]),['index','time','label'])

    H5File("%s/%s.h5"%(outdir,exFn)).toSparseAnnData("%s/%s.h5ad"%(outdir,exFn),BLOCK=5000)
    # # Load in h5 file and convert it to anndata
    d1=anndata.read_h5ad("%s/%s.h5ad"%(outdir,exFn))


    sc.settings.figdir = '%s/figures'%(outdir)
    # # Pre-processing ...
    print("pre-processing...")
    sc.pp.filter_cells(d1,min_genes=200)
    sc.pp.filter_genes(d1,min_cells=3)

    if iformat=='raw':
        MTFlag1=d1.var_names.str.upper().str.startswith('MT-')
        MTFlag2=d1.var_names.str.upper().str.startswith('MT.')
        MTFlag=[bool(a+b) for a,b in zip(MTFlag1,MTFlag2)]
        d1.var['mt'] = MTFlag
        # # plot n_genes, total_counts, and mt counts
        sc.pp.calculate_qc_metrics(d1, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        #sc.pl.violin(d1, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True, show=False, save="_qc.pdf")
        sc.pl.scatter(d1, x='total_counts', y='pct_counts_mt',show=False, save="_mt.pdf")
        sc.pl.scatter(d1, x='total_counts', y='n_genes_by_counts',show=False, save="_n_genes.pdf")
        d1 = d1[d1.obs.pct_counts_mt < 40, :]
        sc.pp.normalize_total(d1, target_sum=1e4)
        sc.pp.log1p(d1)
            
    # # filtering genes based on dispersion
    sc.pp.highly_variable_genes(d1, min_mean=0.0125, max_mean=5, min_disp=mindisp)
    sc.pl.highly_variable_genes(d1,show=False, save=".pdf")
    d1 = d1[:, d1.var.highly_variable]

    # # Removing batch effects
    #sc.pp.regress_out(d1, ['total_counts', 'pct_counts_mt'])
    #sc.pp.scale(d1, max_value=10)

    # # Dimension reduction
    sc.tl.pca(d1, svd_solver='arpack')


    # # Computing the neighborhood graph
    sc.pp.neighbors(d1, n_neighbors=15, n_pcs=50)
    sc.tl.diffmap(d1)

    # # clustering... 
    sc.tl.leiden(d1,resolution=cluRes)
    sc.tl.paga(d1)
    sc.pl.paga(d1,show=False,save="_Traj.pdf")
    sc.tl.umap(d1,init_pos='paga')
    sc.pl.umap(d1,color=['leiden','time'],legend_loc='on data',show=False,save="_clustering.pdf")


    # # get DE genes for each of the clusters
    sc.tl.rank_genes_groups(d1, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(d1, n_genes=25, sharey=False,show=False, save="_global_DE_genes.pdf")


    # # 
    d1.write_h5ad("%s/%s.h5ad"%(outdir,exFn),compression=9)
    print("\n\n>>>>------------------------------------------------<<<<")
    print("prerun completed! please run scdiff2 for the second pass")
    return d1

def main():
    parser=argparse.ArgumentParser(description="scdiff2 pre-run")
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i','--input',required=True,help='input single cell RNA-seq expression data')
    required.add_argument('-o','--output',required=True,help='output directory')
    optional.add_argument('-f','--format',required=False, default='raw', help='the format of input expression, either raw/norm (raw: raw read counts, norm: normalized expression')
    optional.add_argument('--mindisp',required=False,default=0.15,help='the dispersion cutoff to filter genes (genes with dipsersion < this cutoff will be filtered')
    optional.add_argument('--cluRes',required=False, default=1, help="The resolution parameter for the leiden clustering method")
    args = parser.parse_args()
    exFn=args.input
    outdir=args.output 
    iformat=args.format
    mindisp=float(args.mindisp)
    cluRes=float(args.cluRes)
    prerun(exFn,outdir,iformat,mindisp,cluRes)

if __name__=="__main__":
    main()
