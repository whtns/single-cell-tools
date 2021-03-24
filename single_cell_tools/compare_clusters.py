#!/usr/bin/python

import pandas as pd
import os 
import IPython
import re
import sys

indir=sys.argv[1]+"/"
outdir=sys.argv[2]+"/"


def absoluteFilePaths(directory):
	for dirpath,_,filenames in os.walk(directory):
		for f in filenames:
			yield os.path.abspath(os.path.join(dirpath, f))

cluster_dir = indir
csv_files = [x for x in absoluteFilePaths(cluster_dir) if x.endswith("_expression_values.csv")] 

csvs = [pd.read_csv(x, delim_whitespace=True) for x in csv_files]
#~ IPython.embed()
csvs0 = map(lambda df: df[df.p_val < 0.05], csvs)
#~ csvs0 = map(lambda df: df[df.Z < 3.0], csvs)

itx_csvs = pd.concat(csvs0, axis=1, join="inner")
union_csvs = pd.concat(csvs0, axis=0)

#union_csvs.to_csv("diffex_genes_union.csv")

csv_names = map(os.path.basename, csv_files)

csv_names = [i.replace("_stringtie_expression_values", "") for i in csv_names]
csv_names  = [re.sub("^.*trs", "", i) for i in csv_names]

csv_names = [outdir+i for i in csv_names]


diffex_csvs = dict(zip(csv_names, csvs0))

for key,value in diffex_csvs.iteritems():
	value.to_csv(key) 


