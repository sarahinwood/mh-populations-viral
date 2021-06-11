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
library(dplyr)

###########
# GLOBALS #
###########

mh_gc_table <- snakemake@input[["mh_gc_table"]]
coverage_file <- snakemake@input[["coverage_file"]]

####################
## test normality ##
####################

gc_table <- fread(mh_gc_table)
coverage <- fread(coverage_file)

##gc depth table
depth_table <- coverage[,c(1,7)]
gc_depth <- merge(gc_table, depth_table, by.x="#Name", by.y="#rname")
gc_depth$test_group <- tstrsplit(gc_depth$plot_label, " and", keep=c(1))
gc_depth_stats <- gc_depth[order(-test_group)]

##remove other
gc_depth_stats <- subset(gc_depth_stats, !(test_group=="Other contig"))

##test normality of meandepth
shapiro <- shapiro.test(gc_depth_stats$meandepth)
shapiro_chars <- capture.output(print(shapiro))
writeLines(shapiro_chars, con=file(snakemake@output[["shapiro_res"]]))

##cannot use parametric tests

##########################
## non-parametric tests ##
##########################

summary <- group_by(gc_depth_stats, test_group) %>%
  summarise(
    count = n(),
    mean = mean(meandepth, na.rm = TRUE),
    sd = sd(meandepth, na.rm = TRUE),
    median = median(meandepth, na.rm = TRUE),
    IQR = IQR(meandepth, na.rm = TRUE)
  )
summary_chars <- capture.output(print(summary))
writeLines(summary_chars, con=file(snakemake@output[["summary_stats"]]))

##does meandepth content between any groups differ significantly
kruskal_res <- kruskal.test(meandepth ~ test_group, data = gc_depth_stats)
kruskal_chars <- capture.output(print(kruskal_res))
writeLines(kruskal_chars, con=file(snakemake@output[["kruskal_res"]]))

##which groups differ significantly
pair_wilcox_res <- pairwise.wilcox.test(gc_depth_stats$meandepth, gc_depth_stats$test_group,
                                        p.adjust.method = "BH")
wilcox_chars <- capture.output(print(pair_wilcox_res))
writeLines(wilcox_chars, con=file(snakemake@output[["wilcox_res"]]))

#write log
sessionInfo()