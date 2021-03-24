#!/usr/bin/Rscript 

library("optparse")

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-p", "--pyout"), type="character",
              default = NA,
              help=".pkl file output from pseudotime_interactive.py",
              metavar="pyout")
)

opt = parse_args(OptionParser(option_list=option_list))
if (!is.na(opt$pyout)) {
  mypickle <- opt$pyout
} else {
  stop("Please provide .pkl file for analysis! See script usage (--help)")
}

library(reticulate)
# library(seuratTools)
library(seuratTools, lib.loc = "~/rpkgs/devel_install/")

# test0 <- reticulate::py_load_object("single_cell_tools/for_seurat.pkl")
test0 <- reticulate::py_load_object(mypickle)

data <- py_to_r(test0$expression_table) %>% 
  t()

meta_data <- py_to_r(test0$annotation)

pca <- py_to_r(test0$PC_expression) %>% 
  as.matrix()


# rename all files------------------------------
rename_shl <- function(myvec){
  
  myvec <- stringr::str_replace(myvec, "X", "") %>% 
    stringr::str_pad(width = max(nchar(.)), pad = "0")
    
  
  paste0("shl20170407-", myvec)
} 

colnames(data) <- rename_shl(colnames(data))

rownames(pca) <- rename_shl(rownames(pca))

rownames(meta_data) <- rename_shl(rownames(meta_data))

# ------------------------------

colnames(pca) <- paste0("PC_", colnames(pca))

my_seu <- Seurat::CreateSeuratObject(data, meta.data = meta_data)

test1 <- my_seu %>% 
  # Seurat::NormalizeData() %>% 
  Seurat::ScaleData() %>%
  Seurat::FindVariableFeatures() %>%
  Seurat::RunPCA(npcs = 20) %>% 
  # seurat_reduce_dimensions() %>%
  identity()

DimPlot(test1, reduction = "pca")

test1@reductions$pca@cell.embeddings <- pca


myloom = "~/single_cell_projects/sc_RB_devel/20170407-SHL-FACS-Hs_proj/output/velocyto/20170407-SHL-FACS-Hs_proj.loom"
test2 <- velocyto_assay(test1, loom_path = myloom, reduction = "pca")

# ------------------------------

cell.colors <- tibble::as_tibble(test2$cluster, rownames = "cellid") %>%
  tibble::deframe() %>%
  as.factor()

levels(cell.colors) <- scales::hue_pal()(length(levels(cell.colors)))

plot_path <- "output/sunlee_corr_with_ptime/pc_plots_20170407.pdf"
message(glue::glue("saving plot to {plot_path}"))
pdf(plot_path)
velocyto.R::pca.velocity.plot(test2@misc$vel, nPcs = 6, cell.colors = cell.colors, plot.cols = 1, arrow.scale = 0.5)
dev.off()
