#!/usr/local/bin/Rscript

suppressMessages(library(optparse))


# default data ------------------------------------------------------------


# default_expr_mat = "~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/transcripts.tpm_census_matrix.csv"
# default_annotation = "~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/FACS_20171031_sunlee_sample_sheet.csv"

# default_clusters = "~/single_cell_pipeline/scde_input/diffex_by_trs_clusters_1_4/"

#SHL 20171031
default_expr_mat = "~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/FACS_20171031_sunlee_H_sapiens_census_matrix.rds"
default_annotation = "~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/FACS_20171031_sunlee_H_sapiens_census_matrix_meta.csv"
default_cell_info <- "~/single_cell_tools/FACS_1031_2017_SHL_input_files/cells_sets_4.csv"
default_plot_settings <- "~/single_cell_tools/FACS_1031_2017_SHL_input_files/plot_settings_2d.csv"
default_out = "/home/skevin"

#SHL 20170407
# default_expr_mat = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/sunhye_census_matrix_20170407.rds"
# default_annotation = "~/single_cell_pipeline/scde_input/shl_0407_w_centroids_cell_info.csv"
# default_cell_info <- "~/tmp/New_cells_sets_3_1.csv"
# default_plot_settings <- "~/tmp/New_plot_settings_2d.csv"
# default_out = "/home/skevin"

# input_rda <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/shl_0407_plot_diffex_input.rda"
input_rda <- "~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/shl_1031_plot_diffex_input.rda"

# shl <- mget(ls(pattern = "default"))
# save(shl, file = input_rda)
  
if (file.exists(input_rda)){
  load( input_rda)
  list2env(shl, globalenv())
} else {
  default_expr_mat = default_annotation =  default_cell_info =  default_plot_settings = default_out = NA

}


#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-e", "--expr_mat"), type="character", default=default_expr_mat,
              help="expression matrix after census normalization [default= %default]", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=default_annotation,
              help="metadata about cells in input file [default= %default]", metavar="character"),
  make_option(c("-c", "--cellset"), type="character", default=default_cell_info,
              help="tab delimited cell settings file [default= %default]", metavar="character"),
  make_option(c("-p", "--plot_settings"), type="character", default=default_plot_settings,
              help="tab delimited plot settings file [default= %default]", metavar="character"), 
  make_option(c("-o", "--out"), type="character", default=default_out,
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  stop("Please provide all necessary arguments.", call.=FALSE)
}

# load required libraries -------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(gtools))
suppressMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86


# load required functions -------------------------------------------------

find_remove_cells <- function(plot_settings, annotation){
  # browser()
  test <- readLines(plot_settings)
  
  # if (!grepl('remove', test)){
  #   return(NULL)
  # }
  
  vecs <- list()
  mtnames <- c()
  for (i in test){
    if (!grepl("#", i) & grepl('remove', i)){
      lbline = strsplit(i, "\t")
      d = unlist(lbline)[[2]]
      vecs <- append(vecs, lbline)
      # add treatment prefix
      if (d %in% c("shCtrl", "sh733", "sh737")){
        d <- paste0("treatment_", d)
      }

      mtnames <- append(mtnames, d) 
      
    }
  }
  
 

  
  pfx <- tolower(gsub("_.*", "", mtnames))
  valid_pfx  <- which(pfx %in% tolower(colnames(annotation)))
  
  if (length(valid_pfx) == 0){
    return(NULL)
  }
  
  which(pfx %in% annotation$treatment_group)
  
  pfx <- pfx[valid_pfx]
  
  sfx <- gsub(".*_", "", mtnames[valid_pfx])
  
  remove_cells <- purrr::map2(pfx, sfx , function(x, y) annotation[annotation[tolower(x)] == y,])
  
  
  
  remove_cells <- dplyr::bind_rows(remove_cells)
  
  ind <- apply(remove_cells, 1, function(x) all(is.na(x)))
  remove_cells <- remove_cells[ !ind, ]
  
  remove_cells <- unique(remove_cells[,1])
  
}


match_cols <- function(match_vecs, sv_name){
  # browser()
  out=NULL
  for (i in match_vecs){
    vetor <- i
    vetor <- vetor[vetor != ""]
    key <- data.frame(sample_id=vetor[-1], sv_name=rep(gsub(".*_","", vetor[[1]]), (length(vetor)-1)))  
    out <- rbind(out, key)
  }  
  colnames(out) <- c("sample_id", sv_name) 
  return(out)
}

convert_mt_setting <- function(cell_settings, plot_settings){
  # browser()
  test <- readLines(cell_settings)
  # drop blank lines
  test <- test[grepl(".*", test)]
  vecs <- list()
  mtnames <- c()
  for (i in test){
    if (!grepl("#", i)){
      lbline = strsplit(i, "\t")
      d = unlist(lbline)[[1]]
      vecs <- append(vecs, lbline)
      mtnames <- append(mtnames, d)    
    }
    
  }
  # browser()
  # add treatment group label to conform to formatting of other tags
  treatment_ind <- which(mtnames %in% c("shCtrl", "sh733", "sh737"))
  mtnames[treatment_ind] <- paste0("treatment_", mtnames[treatment_ind])
  
  pfx <- unique(gsub("_.*", "", mtnames[grep("_", mtnames)]))
  pfx <- paste0(pfx, "_")
  test <- list()
  for (i in pfx){
    test1 <- list(which(startsWith(mtnames, i)))
    names(test1) = tolower(gsub("_", "", i))
    test <- append(test, test1)
  }
  
  sub_vecs <- list()
  vec_names <- list()
  for (i in test){

    test_vec <- vecs[i]
    sub_vecs <- append(sub_vecs, list(test_vec))
  }
  
  # names(sub_vecs) <- vec_names[1:length(sub_vecs)]
  names(sub_vecs) <- names(test)
  
  sub_vecs <- sub_vecs[unlist(lapply(sub_vecs, length) != 0)]
  
  param_dfs <- purrr::map2(sub_vecs, names(sub_vecs), match_cols)
  
  
  
  if (is.list(param_dfs) & length(param_dfs) != 0) {
    # browser()
    param_dfs <- param_dfs %>%
      Reduce(function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="sample_id"), .) %>% 
      dplyr::arrange(sample_id)  
    
    dup_cells <- which(duplicated(param_dfs[,1]))
    if (any(dup_cells)){
      print(paste0("cells ", paste(param_dfs$sample_id[dup_cells], collapse = " "), " found duplicated in cell sets! They will be removed from analysis"))
      param_dfs <- param_dfs[-dup_cells,]
    }
    
    rownames(param_dfs) <- param_dfs[,1]
    
  }
  
  remove_cells <- find_remove_cells(plot_settings, param_dfs)
  
  
  return(list("annotation" = param_dfs, "removed_cells" = remove_cells))
  
}

take_input <- function(prompt = question, interactive = FALSE){
  if(interactive){
    param <- readline(prompt=question)
  } else {
    param <- readLines("stdin", n=1)
  }
  return(param)
}

lookup_genes <- function(txname){
  
  txs <- transcripts(edb, filter = TxIdFilter(txname),
                     columns = c("symbol"))
  return(txs$symbol)
  
}

lookup_transcripts <- function(genename){
  
  txs <- transcripts(edb, filter = GenenameFilter(genename),
                     columns = c("tx_id"))
  return(txs$tx_id)
  
}

plot_trx_by_treatment_and_facet <- function(transcript, annotation, facet){
  # browser()
  filt_cm <- census_matrix[rownames(census_matrix) %in% transcript,]
  filt_cm <- tidyr::gather(filt_cm, "sample_id", "counts") %>% 
    inner_join(annotation)
  filt_cm <- filt_cm[!is.na(filt_cm[[facet]]),]
  new_levels <- mixedsort(levels(filt_cm[[facet]]))
  filt_cm[[facet]] <- factor(filt_cm[[facet]], levels = new_levels)
  
  RBKD_filt_cm <- filt_cm %>% 
    dplyr::filter(treatment_group %in% c("sh733", "sh737", "sh842")) %>% 
    dplyr::mutate(treatment_group = "RBKD")
  
  comb_filt_cm <- rbind(RBKD_filt_cm, filt_cm)
  
  if(any(comb_filt_cm$day %in% c("day_15"))){
      comb_filt_cm <- dplyr::mutate(comb_filt_cm, day = dplyr::case_when(
        day == "day_3" | day == "day_5" ~ "day 3-5",
        day == "day_7" | day == "day_9" ~ "day 7-9",
        day == "day_12" | day == "day_15" ~ "day 12-15",
        TRUE ~ as.character(day)
      ))
  }
  
  comb_filt_cm <- comb_filt_cm[!is.na(comb_filt_cm$day),]
  
  comb_filt_cm$day <- factor(comb_filt_cm$day, levels <- unique(comb_filt_cm$day))
  # comb_filt_cm$day <- factor(comb_filt_cm$day, levels = c("RBKD", levels(filt_cm$treatment_group)))
  # comb_filt_cm$treatment_group <- factor(comb_filt_cm$treatment_group, levels = c("RBKD", levels(filt_cm$treatment_group)))

  bplot <- ggplot(data = comb_filt_cm, aes_string(x=facet, y="counts")) + 
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.1) +
    # scale_x_discrete(mixedsort(levels(filt_cm[[facet]]))) +
    facet_grid(. ~ treatment_group) + 
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    labs(title = lookup_genes(transcript), subtitle = transcript)

  return(bplot)
}

plot_genes_summed_trx <- function(census_matrix, transcripts, annotation, facet){
  # browser()
  # input = readline(prompt="Enter gene symbol to plot (ex. MYCN): ")
  # prep_transcripts <- lapply(transcripts, prep_transcript)
  bplots <- lapply(transcripts, plot_trx_by_treatment_and_facet, annotation, facet)
  return(bplots)
}

# load data ---------------------------------------------------------------

print("loading census matrix")
# census_matrix = cataract::safe_read(opt$expr_mat)
census_matrix = readRDS(opt$expr_mat)
print("loading cell metadata")
annotation = cataract::safe_read(opt$annotation)

mt_settings <- convert_mt_setting(opt$cellset, opt$plot_settings)

if (!is.null(mt_settings$annotation)) {
  # removed_print <- annotation[annotation[,"sample_id"] %in% mt_settings$removed_cells,]
  # print("removing cells")
  # print(removed_print)
  
  exist_cols <- colnames(mt_settings$annotation)[!colnames(mt_settings$annotation) %in% colnames(annotation)]
  exist_cols <- c("sample_id", exist_cols)
  
  mt_settings$annotation <- dplyr::select(mt_settings$annotation, exist_cols)
  # mt_settings$annotation <- mt_settings$annotation[complete.cases(mt_settings$annotation),]
  oldannotation <- as.data.frame(annotation)
  names(oldannotation) <- tolower(names(oldannotation))
  annotation <- left_join(oldannotation, mt_settings$annotation, by = "sample_id")
  
  #remove cells from annoation and expression matrix
  annotation <- annotation[!annotation$sample_id %in% mt_settings$removed_cells, ]
  census_matrix <- census_matrix[,colnames(census_matrix) %in% annotation$sample_id] 
}


# run while loop ----------------------------------------------------------


while (TRUE) {
  question = "Choose from following:\n
                                    [P] Plot gene by treatment and day\n
                                    [M] Plot gene by treatment and cluster\n
                                    [L] Plot list of genes by treatment and day\n
                                    [C] Plot list of genes by treatment and cluster\n
                                    [X] Exit \n"
  
  cat(question)
  action <- toupper(take_input(question))
  # action <- readline(prompt=question)
  
  if(action == "X"){
    break
  } else if (action == "P"){
    
    # plot gene by treatment and day ---------------------------------
    
    gene_question <- "Enter gene symbol to plot (ex. MYCN): "
    cat(gene_question)
    input <- take_input(gene_question)
    transcripts <- lookup_transcripts(input)
    transcripts <- transcripts[which(transcripts %in% rownames(census_matrix))]
    bplots <- plot_genes_summed_trx(census_matrix, transcripts, annotation, "day") 
    pdf_out <- paste0(opt$out, "/", input, ".pdf")
    pdf(pdf_out)
    invisible(lapply(bplots, print))
    dev.off()
    print(paste0("saving ", pdf_out))
  } else if (action == "M"){
    
    # plot gene by treatment and cluster ---------------------------------
    
    gene_question <- "Enter gene symbol to plot (ex. MYCN): "
    cat(gene_question)
    input <- take_input(gene_question)
    transcripts <- lookup_transcripts(input)
    transcripts <- transcripts[which(transcripts %in% rownames(census_matrix))]
    bplots <- plot_genes_summed_trx(census_matrix, transcripts, annotation, "cluster") 
    pdf_out <- paste0(opt$out, "/", input, ".pdf")
    pdf(pdf_out)
    invisible(lapply(bplots, print))
    dev.off()
    print(paste0("saving ", pdf_out))
  } else if (action == "L"){
    

  # Plot list of genes by treatment and day -----------------------------

    treatment_question <- paste0("Provide filepath to list of genes (as .txt file)")
    
    cat(treatment_question)
    gene_list <- take_input(treatment_question)
    # diffex_group <- readline(prompt=treatment_question)
    gene_list <- read.csv(gene_list)
    
    # top_n_question <- "how many genes do you want to plot? "
    # cat(top_n_question)
    # top_n <- take_input(top_n_question)
    
    # transcripts <- rownames(diffex_genes[c(1:top_n),])
    
    transcripts <- gene_list[,1]
    
    transcripts <- transcripts[which(transcripts %in% rownames(census_matrix))]
    bplots <- plot_genes_summed_trx(census_matrix, transcripts, annotation, "day") 
    pdf_out <- paste0(opt$out, "/", diffex_group, ".pdf")
    pdf(pdf_out)
    invisible(lapply(bplots, print))
    dev.off()
    print(paste0("saving ", pdf_out))
  } else if (action == "C"){
    

  # Plot list of genes by treatment and cluster -----------------------------
    
    treatment_question <- paste0("Provide filepath to list of genes (as .txt file)")
    
    cat(treatment_question)
    diffex_group <- take_input(treatment_question)
    # diffex_group <- readline(prompt=treatment_question)
    gene_list <- read.csv(gene_list)
    
    # top_n_question <- "how many genes do you want to plot? "
    # cat(top_n_question)
    # top_n <- take_input(top_n_question)
    
    transcripts <- rownames(diffex_genes[c(1:top_n),])
    
    transcripts <- transcripts[which(transcripts %in% rownames(census_matrix))]
    bplots <- plot_genes_summed_trx(census_matrix, transcripts, annotation, "cluster") 
    pdf_out <- paste0(opt$out, "/", basename(diffex_group), ".pdf")
      pdf(pdf_out)
    invisible(lapply(bplots, print))
    dev.off()
    print(paste0("saving ", pdf_out))
  } else if (action == "I"){
    

  # debug -------------------------------------------------------------------

    browser()
    print("a")
  } 
}

