#!/usr/bin/env python
"""
author: jun ding
date: 2020-07-08

function: this code is used to generate  a tab-formated expression file required by scdiff software suite 

"""
import sys,os,pdb
import argparse 
import scipy.io
import csv
from scipy.sparse import csr_matrix
import gzip

def main():
	parser=argparse.ArgumentParser(description="convert 10x genomics mtx.matrix, barcodes.csv, and genes.tsv to a tab expression file that scdiff2 accepts")
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')
	required.add_argument('-i','--inputFolder',required=True,help='the input folder that stores mtx.matrix, barcodes.csv, and genes.tsv' +
																'you can also provide an optional meta.tsv to provide the time and label information for each of cells' +
																'if not provided, by default, all the cells will be set as time point 1 and label NA')
	args = parser.parse_args()
	exFn=args.inputFolder.strip()
	exFn=exFn.split("/")
	exFn=[item for item in exFn if item!=""]
	exFn="/".join(exFn)
	try:
		mat=scipy.io.mmread(f"{exFn}/matrix.mtx.gz")
		mat=mat.transpose()
		feature_ids = [row[0] for row in csv.reader(gzip.open(f"{exFn}/genes.tsv.gz","rt"), delimiter="\t")]
		gene_names = [row[1] for row in csv.reader(gzip.open(f"{exFn}/genes.tsv.gz","rt"), delimiter="\t")]
		barcodes = [row[0] for row in csv.reader(gzip.open(f"{exFn}/barcodes.tsv.gz","rt"), delimiter="\t")]
	except:
		mat=scipy.io.mmread(f"{exFn}/matrix.mtx")
		mat=mat.transpose()
		gene_names = [row[1] for row in csv.reader(open(f"{exFn}/genes.tsv","rt"), delimiter="\t")]
		barcodes = [row[0] for row in csv.reader(open(f"{exFn}/barcodes.tsv","rt"), delimiter="\t")]
	
	mat=csr_matrix(mat)
	if os.path.exists(f"{exFn}/meta.tsv.gz"):
		meta=[row[1] for row in csv.read(gzip.open(f"{exFn}/meta.tsv.gz","rt"),delimiter="\t")]
	elif os.path.exists(f"{exFn}/meta.tsv"):
		meta=[row for row in csv.read(open(f"{exFn}/meta.tsv","rt"),delimiter="\t")]
	else:
		meta=[]
	
	FR=['cell','time','label']+gene_names
	#pdb.set_trace()
	g=open(exFn+".E","a")
	g.write("\t".join(FR)+"\n")
	
	for i in range(len(barcodes)):
		ci=barcodes[i]
		ti=meta[i][1] if len(meta)>0 else 1
		li=meta[i][2] if len(meta)>0 else 'NA'
		ei=mat[i].toarray()[0]
		oi=[ci,ti,li]+list(ei)
		oi=[str(item) for item in oi]
		g.write("\t".join(oi)+"\n")
		print(i)
	g.close()
	
if __name__=='__main__':
	main()
	
