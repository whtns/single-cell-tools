#!/usr/bin/env python 

import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np
import scipy as sc
import sklearn
from sklearn import decomposition
from sklearn import cluster
import seaborn as sns
import math
import IPython
sys.path.insert(0, "/home/hp/CHLA/single-cell-tools/")
from sc_pseudotime import *
from matplotlib.backends.backend_pdf import PdfPages
import os
import ipdb 
import pickle
import functools

from inspect import currentframe, getframeinfo

import argparse

"~/python_packages/single-cell-tools/resources/input/FACS_0407_2017_SHL_input_files/sunhye_census_matrix_20170407-comma_delimited.csv"

parser = argparse.ArgumentParser(description="runs genes_correlated_with_pseudotime")
parser.add_argument("-e", "--expression-matrix", dest="expr_mat", default="/home/skevin/single_cell_projects/sc_RB_devel/20170407-SHL-FACS-Hs_proj/results/sunhye_census_matrix_20170407.csv", help="gene by cell matrix of expression values", metavar="EXPR")
parser.add_argument("-c", "--cell-sets", dest="cell_sets", default="/home/skevin/python_packages/single-cell-tools/resources/sunlee_input/New_cells_sets_3_6.csv", help="cell sets", metavar="CELL_SETS")
parser.add_argument("-p", "--plot-settings", dest="plot_settings", default="/home/skevin/python_packages/single-cell-tools/resources/sunlee_input/New_plot_settings_2d.csv", help="plot settings", metavar="PLOT_SETTINGS")
parser.add_argument("-r", "--corr-method", dest="corr_method", default="spearman", help="method of correlation (spearman or pearson)", metavar="CORR_METHOD")
parser.add_argument("-f", "--feature", dest="feature", default = "g", help="feature of interest; either 'gene' or 'transcript' depending on desired output", metavar="FEATURE", required=False)
parser.add_argument("-o", "--outdir", dest="outdir", default="../output/sunlee_corr_with_ptime", help="a name to give to the output file", metavar="OUTDIR")
parser.add_argument("-pt", "--pseudotime", dest="pseudotime", default = "../resources/example_input_files/PT1_shCtrl.csv", help="experimental cells. a list of pseudotime values. Can accept multiple values", metavar="PTIME", nargs='+', required=False)
parser.add_argument("-cpt", "--control-pseudotime", dest="ctrl_pseudotime", help="control cells. a list of pseudotime values. Can accept multiple values", metavar="PTIME", nargs='+') # default = "../resources/example_input_files/shCtrl_PT.csv"

# default = "resources/example_input_files/shCtrl_PT.csv"

try:
    options = parser.parse_args()
except SystemExit as err: 
	if err.code == 2: 
		parser.print_help()
		sys.exit(0)
 
expression_file = os.path.expanduser(options.expr_mat)
cellset_file    = os.path.expanduser(options.cell_sets)
settings_file   = os.path.expanduser(options.plot_settings)
correlation_method = options.corr_method
output_dir      = options.outdir+"/"

if not os.path.exists(output_dir):
  print(output_dir+"doesn't exist yet. Creating "+output_dir+" first")
  os.mkdir(output_dir)

# pseudotime_files = []
# pseudotime_files.append(options.pseudotime)
pseudotime_files = sorted(options.pseudotime)
print(pseudotime_files)

if options.ctrl_pseudotime:
  # ctrl_pseudotime_files = []
  # ctrl_pseudotime_files.append(options.ctrl_pseudotime)
	ctrl_pseudotime_files = sorted(options.ctrl_pseudotime)
sett = settings(settings_file, cellset_file)
expression_table, annotation = read_expression(expression_file, sett, min_expression = 0.1, min_cells = 5)
#~ expression_table = trx_expression_table.copy()

def symbols_from_geneids(geneids):
	mg = mygene.MyGeneInfo()
	# ~ for i in enumerate(genes): 
	gene_info = mg.querymany(geneids.index, scopes='ensembl.gene', fields='symbol')
	gene_info[:] = [d for d in gene_info if d.get('notfound') != True]
	symbols = [d['symbol'] for d in gene_info]
	return(symbols)

gene_expression_file = output_dir+"gene_expression.csv"
if not os.path.isfile(gene_expression_file):
	print("creating gene expression table; saving as "+gene_expression_file)
	gene_trx_dic = get_gene_transcript_dic(expression_table)
	gene_expression_table =  trx_to_gene_exp_table(expression_table, gene_trx_dic)
	gene_expression_table.to_csv(gene_expression_file, sep = "\t")
	save_obj(output_dir, gene_trx_dic, 'gene_trx_dic')
else:
	gene_expression_table = pd.read_csv(gene_expression_file, index_col=0, sep ="\t")
	gene_trx_dic = load_obj(output_dir, 'gene_trx_dic')
	
# 

# @functools.lru_cache(maxsize=None)
def set_pts(pseudotime_files, cell_set_flag):
  pt = [read_pseudotime_from_file(i) for i in pseudotime_files]
  corr = [get_correlation_with_pseudotime(x, expression_table, annotation, gene_trx_dic, cell_set_flag, feature = options.feature, method=correlation_method) for x in pt]
  #~ corr = [corr for i,corr in enumerate(d['spearman'] for d in corr_exp_dict)]
  corr = pd.concat(corr, axis=1)
  ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in pseudotime_files]
  pt = dict(zip(ptime_titles, pt))
  #~ corr_pt_dict = {'corr': corr, 'pt': pt, 'exp': corr_exp_dict[0]['exp']}
  return corr, pt
  


#load ctrl pseudotime objects if control files supplied
try:
	ctrl_pseudotime_files
except:
	print("no ctrl pseudotime files supplied!")
	ctrl_user_ptimes = "none"
	cpt = "none"
else:
	# read in control pseudotime files
	cpt = map(read_pseudotime_from_file, ctrl_pseudotime_files)
	cptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in ctrl_pseudotime_files]
	cpt = dict(zip(cptime_titles, cpt))
	if (options.feature == "t" or options.feature == "transcript"):
		ctrl_correlation_file = "_".join(cptime_titles)+"_"+correlation_method+"_transcript_correlation.csv"
	elif (options.feature == "g" or options.feature == "gene"):
		ctrl_correlation_file = "_".join(cptime_titles)+"_"+correlation_method+"_symbol_correlation.csv"

	#~ gene_exp_file = "_".join(cptime_titles)+"_"+correlation_method+"_gene_expression.csv"
	
	# read correlation files from similarly named files
	if os.path.exists(ctrl_correlation_file):
		print(ctrl_correlation_file)
		ctrl_corr = pd.read_csv(ctrl_correlation_file, sep="\t", index_col=0)
	
	# check if control correlation files have already been read in
	try:
		ctrl_corr
	except:
		ctrl_corr, pt = set_pts(ctrl_pseudotime_files, cell_set_flag="ctrl")
		ctrl_corr.columns = sorted(cpt.keys())
		#~ exp.to_csv(gene_exp_file, sep="\t")
		ctrl_corr.to_csv(ctrl_correlation_file, sep="\t")
	else:
		if sorted(cpt.keys()) == list(ctrl_corr.columns):
			pass
		else:
			print("column names do not match!")
			
	ctrl_user_ptimes = ' '.join(cpt.keys())


#load experimental pseudotime objects (required)

# read in pseudotime files
pt = map(read_pseudotime_from_file, pseudotime_files)
ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in pseudotime_files]
pt = dict(zip(ptime_titles, pt))

if options.feature == "t":
	correlation_file = "_".join(ptime_titles)+"_"+correlation_method+"_transcript_correlation.csv"
elif options.feature == "g":
	correlation_file = "_".join(ptime_titles)+"_"+correlation_method+"_symbol_correlation.csv"

gene_exp_file = "_".join(ptime_titles)+"_"+correlation_method+"_gene_expression.csv"

print(correlation_file)


if os.path.isfile(gene_exp_file):
	exp = pd.read_csv(gene_exp_file, sep = "\t", index_col=0)


# read correlation files from similarly named files
if os.path.exists(correlation_file):
	corr = pd.read_csv(correlation_file, sep="\t", index_col=0)

corr_columns = []
for i in sorted(pt.keys()):
	corr_columns += [i+"_exp_corr"]
	corr_columns += [i+"_ctrl_corr"]

# check if experimental correlation files have already been read in
try:
	corr
except:
	if cpt == "none":
		corr, pt = set_pts(pseudotime_files, cell_set_flag="exp")
		corr_columns = []
		for i in sorted(pt.keys()):
			corr_columns += [i+"_exp_corr"]
		corr.columns = corr_columns
		corr.to_csv(correlation_file, sep="\t")
	else:
		corr, pt = set_pts(pseudotime_files, cell_set_flag="mix")
		corr_columns = []
		for i in sorted(pt.keys()):
			corr_columns += [i+"_exp_corr"]
			corr_columns += [i+"_ctrl_corr"]
		corr.columns = corr_columns
		corr.to_csv(correlation_file, sep="\t")
else:
	if corr_columns == list(corr.columns):
		pass
	else:
		print("column names do not match!")


user_ptimes = ' '.join(pt.keys())

## function finds genes within threshold
def genes_within_threshold(corr, ht, lt=None):
  if lt is not None:
    corr = corr[corr.ctrl < lt]
    corr = corr[corr.exp > ht]
    corr = corr[['exp', 'ctrl']]
    corr = corr.round(4)
    corr = corr.reset_index(col_fill="gene")
    corr = corr.sort_values(by=['exp', 'ctrl'], ascending=[False, True])
  else:
    corr = corr[corr.exp > ht]
    corr = corr[['exp']]
    corr = corr.round(4)
    corr = corr.reset_index(col_fill="gene")
    corr = corr.sort_values(by='exp', ascending=False)
  
  genes_of_interest = list(zip(*map(corr.get, corr)))

  
  return(genes_of_interest)

## function plots genes of interest (pd.Index) into pdf
def plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, ptime, pt, cpt=None, squeeze=True):
  # IPython.embed()
  
  pt_ids = {'pt' : list(pt.keys())}
  if cpt is not None:
    pt.update(cpt)
    pt_ids['cpt'] = list(cpt.keys())
  
  # IPython.embed()
  gene_plots = [plot_gene_with_pseudotime(expression_table, pt, pt_ids, gene_tuple, annotation) for gene_tuple in genes_of_interest]
    
  # IPython.embed()
  save_as_pdf_pages(gene_plots, out_filename)
  

#~ default_output_dir = "genes_corr_w_ptime"
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

if options.feature == "g":
	expression_table = gene_expression_table


try:
	ctrl_corr
except:
  total_corr = corr
else:
  total_corr = corr.append(ctrl_corr)

# main loop / choosing action
while True:
  question = """Choose from following:
  [C]	Create Correlation file from New Pseudotime Files
  [D]	Plot User-Supplied Genes
  [T]	Plot Top N Features (Genes or Transcripts) with highest correlation
  [X]	Exit
  """
  action = input(question).upper()
  if(action == "X"):
  	break
  	#~ exit()
  	#~ FACS_0407_2017_SHL_input_files/DEGS_day_12.csv
  elif(action == "C"):
  	ptime_paths = input("provide path(s) to new pseudotime files ").split(",")
  	corr_out = input("provide filename for new correlation files ")
  	## block of code to calculate correlations
  	pt = map(read_pseudotime_from_file, ptime_paths)
  	corr = [get_correlation_with_pseudotime(x, expression_table, gene_trx_dic, method=correlation_method) for x in pt]
  	corr.to_csv(corr_out+"_"+correlation_method+"_correlation.csv", sep="\t")
  	ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in ptime_paths]
  	ptimes = dict(zip(ptime_titles, ptime_paths))
  	user_ptimes = ' '.join(ptime_titles)
  	corr.columns = ptime_titles
  	
  elif(action == "D"):
  	DEG_path = input("provide path to differentially expressed genes ")
  	ptime = input("Which pseudotime would you like correlate with? ("+user_ptimes+ ") ")
  	ctrl_ptime = input("Which ctrl pseudotime would you like to correlate with? ("+ctrl_user_ptimes+ ") ")
  	DEGS = pd.read_csv(DEG_path, index_col=0, header=None)
  	if options.feature == "g":
  		DEGS = symbols_from_geneids(DEGS)
  	elif options.feature == "t":
  		DEGS = DEGS.index
  	corr["order"] = corr[ptime+"_exp_corr"].abs()
  
  	DEGS = corr[corr.index.isin(DEGS)].index
  	out_filename = output_dir+correlation_method+"_"+ptime+"_DEGS.pdf"
  	# 
  	if ctrl_ptime == '':
  		if len(pt) == 1:
  		  # 
  		  plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, ptime, pt, squeeze=False)
  			
  		else:
  			plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, ptime, pt)

  	else:
  		plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, ptime, pt, cpt)
  		# plot transcripts of interest
  		#~ plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, ptime, pt, cpt[ctrl_ptime])
  elif(action == "T"):
    top_n = int(input("How many genes would you like to plot? "))
    ptime = input("Which pseudotime would you like to order by? ("+user_ptimes+ ") ")
    ctrl_ptime = input("Which ctrl pseudotime would you like to correlate with? Leave blank if no control present. ("+ctrl_user_ptimes+ ") ")

    ht = input("set upper threshold (def. 0.3) ")
    if not ht == '':
      ht = float(ht)
    else:
      ht = 0.3
    lt = float(input("set lower threshold (def. 0.2) Leave blank if no control present. ") or 0)
    corr["exp"] = (corr[ptime+"_exp_corr"]).abs()
    
    # IPython.embed()
    if ctrl_ptime != '':
      corr["ctrl"] = ctrl_corr[ctrl_ptime].abs()
      genes_of_interest = genes_within_threshold(corr, ht, lt)
    else:
      genes_of_interest = genes_within_threshold(corr, ht)
    
    if top_n > len(genes_of_interest):
    	print("error! number of genes requested exceeds number of genes matching filtering criteria ("+str(len(genes_of_interest))+")")
    	pass
    	
    if ctrl_ptime == '':
    	genes_of_interest = genes_of_interest[0:top_n]
    	out_filename = output_dir+correlation_method+"_"+ptime+"_top_"+str(top_n)+"_genes.pdf"
    	plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, ptime, pt)
    else:
      # IPython.embed()
      genes_of_interest = genes_of_interest[0:top_n]
      out_filename = output_dir+correlation_method+"_"+ptime+"_top_"+str(top_n)+"_genes.pdf"
      plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, ptime, pt, cpt)
  elif(action == "I"):
    IPython.embed()

#~ if __name__ == "__main__":
	#~ main()
