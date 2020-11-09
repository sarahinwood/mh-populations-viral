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

st_depth_file <- snakemake@input[["samtools_depth"]]
##scaffold labels
scaffold_id_table <- snakemake@input[["scaffold_id_table"]]

########
# MAIN #
########

st_depth_names <- c("Scaffold_full_id", "BP", "depth")
st_depth <- fread(st_depth_file, col.names=st_depth_names)
##merge with table of scaffold ids for plotting
scaffold_table <- fread(scaffold_id_table, header=TRUE)
st_depth_boxpl <- merge(st_depth, scaffold_table, by="Scaffold_full_id", all.x=TRUE)
#order chr numerically
st_depth_boxpl$Scaffold_id <- factor(st_depth_boxpl$Scaffold_id,
    levels=c("PGA_scaffold0", "PGA_scaffold1", "PGA_scaffold2",
      "PGA_scaffold3", "PGA_scaffold4", "PGA_scaffold5",
      "PGA_scaffold6", "PGA_scaffold7", "PGA_scaffold8",
      "PGA_scaffold9", "PGA_scaffold10", "PGA_scaffold11",
      "Scaffold1791", "Scaffold5674", "Scaffold6852",
      "Scaffold13392", "Scaffold15344", "Scaffold15995",
      "Scaffold16687", "Scaffold18966", "Scaffold27879",
      "Scaffold28498", "Scaffold29164", "Scaffold30641", "Scaffold32600", 
      ##non-LbFV containing contigs
      "Scaffold3332", "Scaffold8243", "Scaffold9703",
      "Scaffold32814", "Scaffold43528", "Scaffold45974"))
##order legend labels
st_depth_boxpl$Scaffold_label <- factor(st_depth_boxpl$Scaffold_label,
  levels=c("Hi-C Scaffold", "Viral Contig, LbFV hit", "Viral Contig"))

###################
##making box plot##
###################

pdf(snakemake@output[["boxplot_y_zoom"]])
ggplot(st_depth_boxpl, aes(x=Scaffold_id, y=depth, colour=Scaffold_label))+
  ##make outlier points somewhat transparent
  geom_boxplot(outlier.shape=NA)+
  theme_light(base_size=18)+
  ##turn scaffold labels 90 degrees, align with middle of tick up against it, fill facte titles darker grey
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun.y=mean, geom="point", colour="grey35")+
  ##colour-blind friendly palette
  scale_colour_viridis(discrete=TRUE)+
  ## - full plot y axis range so large boxplots not visible, just outliers
  coord_cartesian(ylim = c(0, 30))
dev.off()

pdf(snakemake@output[["boxplot"]])
ggplot(st_depth_boxpl, aes(x=Scaffold_id, y=depth, colour=Scaffold_label))+
  ##make outlier points somewhat transparent
  geom_boxplot()+
  theme_light(base_size=18)+
  ##turn scaffold labels 90 degrees, align with middle of tick up against it, fill facte titles darker grey
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun.y=mean, geom="point", colour="grey35")+
  ##colour-blind friendly palette
  scale_colour_viridis(discrete=TRUE)
dev.off()

#write log
sessionInfo()
