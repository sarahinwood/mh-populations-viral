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
library(ggplot2)
library(viridis)

###########
# GLOBALS #
###########
##depth files
st_depth_file_i12 <- snakemake@input[["st_depth_file_i12"]]
st_depth_file_i28 <- snakemake@input[["st_depth_file_i28"]]
st_depth_file_i52 <- snakemake@input[["st_depth_file_i52"]]
st_depth_file_i84 <- snakemake@input[["st_depth_file_i84"]]
##scaffold labels
mh_gc_table <- snakemake@input[["mh_gc_table"]]

########
# MAIN #
########

## read in depth tables ##
st_depth_names <- c("#Name", "BP", "depth")
st_depth_12 <- fread(st_depth_file_i12, col.names=st_depth_names)
st_depth_28 <- fread(st_depth_file_i28, col.names=st_depth_names)
st_depth_52 <- fread(st_depth_file_i52, col.names=st_depth_names)
st_depth_84 <- fread(st_depth_file_i84, col.names=st_depth_names)
##add sample label so tables can be joined
st_depth_12$Sample <- "Pupa 1"
st_depth_28$Sample <- "Pupa 2"
st_depth_52$Sample <- "Larva"
st_depth_84$Sample <- "Pupa 3"
full_depth_table <- rbind(st_depth_12, st_depth_28, st_depth_52, st_depth_84)

## scaffold ID table ##
scaffold_table <- fread(mh_gc_table, header=TRUE)
scaffold_table$plot_group <- tstrsplit(scaffold_table$plot_label, " and", keep=c(1))
scaffold_table$scaffold_number <- tstrsplit(scaffold_table$'#Name', "_", keep=c(2))
scaffold_table$scaffold_number <- as.numeric(as.character(scaffold_table$scaffold_number))
setorder(scaffold_table, scaffold_number)

## full table for plotting ##
st_depth_labels <- merge(full_depth_table, scaffold_table, by="#Name", all.x=TRUE)
st_depth_boxpl <- subset(st_depth_labels, !(plot_group=="Other contig"))
##remove outlier contigs
mh_depth_outliers <- list("scaffold_90", "scaffold_995")
st_depth_boxpl <- subset(st_depth_labels, !(`#Name` %in% mh_depth_outliers))
##order legend labels
st_depth_boxpl$plot_group <- factor(st_depth_boxpl$plot_group, levels=c("Hi-C scaffold", "Viral contig"))

##########################
##making panel box plot##
##########################
##height/width in inches
pdf(snakemake@output[["boxplot_panel"]], height=7.5, width=10)
ggplot(st_depth_boxpl, aes(x=reorder(`#Name`, scaffold_number), y=depth, colour=plot_group))+
  geom_boxplot(outlier.shape=NA)+
  theme_bw(base_size=18)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank())+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun.y=mean, geom="point", colour="grey35")+
  scale_colour_viridis(discrete=TRUE)+
  coord_cartesian(ylim = c(0,35))+
  facet_wrap(~Sample)
dev.off()

#write log
sessionInfo()
