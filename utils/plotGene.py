#!/usr/bin/env python
"""
author: Jun Ding
date: 2020-07-06

function: plot the expression of input gene

copy and modification of this code is allowed for academic purposes.  
Please don NOT remove this author statement under any condition.  
"""

import sys,os,pdb,argparse
import anndata 
import scanpy as sc

def plotGene(exFn,gene):
	prRes=anndata.read_h5ad(exFn)
	if gene in prRes.var.index:
		sc.pl.umap(prRes,color=[gene])
	else:
		print("Error! please check your input gene ID, it must be the same as in your expression file")  
		print("Also, the missing gene could be caused by the dispersion based gene filtering by the prerun program")
	
def main():
	parser=argparse.ArgumentParser(description="scdiff2 plotGene")
	parser.add_argument('-i','--input',required=True,help='input h5ad prerun result')
	parser.add_argument('-g','--gene',required=True, help='gene name you want to explore, must be the same ID as in your original input expression file')

	args = parser.parse_args()
	exFn=args.input
	gene=args.gene 
	plotGene(exFn,gene)


if __name__=="__main__":
	main()
