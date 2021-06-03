#!/usr/bin/Rscript

library(gtools)

# Read in the data
x <- scan("~/single_cell_tools/FACS_0407_2017_SHL_input_files/cluster info.csv", what="", sep="\n")
# Separate elements by one or more whitepace
y <- strsplit(x, "\t")
y <- lapply(y, function(a) a[a!=""])
# Extract the first vector element and set it as the list element name
names(y) <- sapply(y, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
y <- lapply(y, `[`, -1)
#y <- lapply(y, function(x) x[-1]) # same as above

df_names = unique(gsub("_.*","", names(y)))


grep(df_names[1], names(y))

vec_sets = lapply(df_names, grep, names(y))

dfs = lapply(vec_sets, function(x) stack(y[x]))

df <- dfs %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="values"), .) %>% 
  transmute(values = values,
            centroid_737 = gsub(".*_", "", ind.x),
            centroid_ctrl = gsub(".*_", "", ind.y),
            centroid_733_234 = gsub(".*_", "", ind.x.x),
            centroid_733_345 = gsub(".*_", "", ind.y.y))

df <- df[mixedorder(df$values),]

cell_info_path = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/scde_input/shl_0407_cell_info.csv"
new_path = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/scde_input/shl_0407_w_centroids_cell_info.csv"

test = read.table(cell_info_path, sep = "\t", header = TRUE)

new_test = left_join(test, df, by = c("Sample_ID" = "values"))

write.table(new_test, new_path, sep = "\t")

