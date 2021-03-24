#!/usr/bin/env python
# coding: utf-8

import loompy
import scvelo as scv
import anndata

with open ("for_seurat.pkl", "rb") as f:
    script_pkl = pickle.load(f)

def rename_shl(mylist):
    mylist = [i.replace("X", "") for i in mylist]
    max_nchar = len(max(mylist, key = len))
    mylist = [i.zfill(max_nchar) for i in mylist]
    mylist = ["shl20170407-"+i for i in mylist]
    return mylist

def plot_velocity(script_pkl, loom_path, components='1,2'):

    data = script_pkl["expression_table"]
    pca = script_pkl["PC_expression"]
    meta_data = script_pkl["annotation"]

    data.index = rename_shl(data.index)
    pca.index = rename_shl(pca.index)
    meta_data.index = rename_shl(meta_data.index)

    adata = anndata.AnnData(data, meta_data)

    # Crate an analysis object
    adata_loom = scv.read("../20170407-SHL-FACS-Hs_proj.loom", cache = True)

    retained_cells = list(set(adata_loom.obs.index).intersection(set(adata.obs.index)))
    retained_cells.sort()
    adata_loom = adata_loom[retained_cells,:]

    adata_loom.var_names_make_unique()

    # plot proporations spliced/unspliced
    # scv.pl.proportions(adata_loom)


    # #preprocess

    scv.pp.filter_and_normalize(adata_loom, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata_loom, n_pcs=30, n_neighbors=30)

    adata_loom.obsm['X_pca'][:,0:20] = pca

    scv.tl.velocity(adata_loom)
    scv.tl.velocity_graph(adata_loom)

    scv.pl.velocity_embedding(adata_loom, basis='pca', components=components)
