#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)

###########
# GLOBALS #
###########

depth <- snakemake@input[["depth"]]

########
# MAIN #
########

depth_table <- fread(depth, header=FALSE)
scaffold_ids <- unique(depth_table$`V1`)

MEAN_SCAFFOLD_DEPTH <- function(x, depth){
  scaffold_depth<-depth[depth$`V1`==x, mean(V3)]
  return(data.table(scaffold_id=x, mean_depth=scaffold_depth))
}

scaffold_depths <- lapply(scaffold_ids, MEAN_SCAFFOLD_DEPTH, depth=depth_table)
scaff_depths_table <- rbindlist(scaffold_depths)
fwrite(scaff_depths_table, snakemake@output[["mean_depth_table"]])

#write log
sessionInfo()
