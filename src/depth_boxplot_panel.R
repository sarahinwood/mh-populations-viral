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
scaffold_id_table <- snakemake@input[["scaffold_id_table"]]

########
# MAIN #
########

##read in depth tables
st_depth_names <- c("Scaffold_full_id", "BP", "depth")
st_depth_12 <- fread(st_depth_file_i12, col.names=st_depth_names)
st_depth_28 <- fread(st_depth_file_i28, col.names=st_depth_names)
st_depth_52 <- fread(st_depth_file_i52, col.names=st_depth_names)
st_depth_84 <- fread(st_depth_file_i84, col.names=st_depth_names)
##add sample label so tables can be joined
st_depth_12$Sample <- "i12_Lincoln"
st_depth_28$Sample <- "i28_Lincoln"
st_depth_52$Sample <- "i52_Lincoln"
st_depth_84$Sample <- "i84_Lincoln"
full_depth_table <- rbind(st_depth_12, st_depth_28, st_depth_52, st_depth_84)
##merge with table of scaffold ids for plotting
scaffold_table <- fread(scaffold_id_table, header=TRUE)
st_depth_boxpl <- merge(full_depth_table, scaffold_table, by="Scaffold_full_id", all.x=TRUE)
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

##########################
##making panel box plot##
##########################
##height/width in inches
pdf(snakemake@output[["boxplot_panel"]], height=15, width=18)
ggplot(st_depth_boxpl, aes(x=Scaffold_id, y=depth, colour=Scaffold_label))+
  ##make outlier points somewhat transparent
  geom_boxplot(outlier.shape=NA)+
  theme_light(base_size=18)+
  ##turn scaffold labels 90 degrees, align with middle of tick up against it, fill facte titles darker grey
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    strip.background=element_rect(fill="grey24"))+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun.y=mean, geom="point", colour="grey35")+
  ##colour-blind friendly palette
  scale_colour_viridis(discrete=TRUE, direction=-1)+
  ##make plot with panel for each sample
  facet_wrap(~Sample)+
  ## - full plot y axis range so large boxplots not visible, just outliers
  coord_cartesian(ylim = c(0, 30))
dev.off()

##trial pointsize sometime to improve size of text etc
#pdf("15_22_pointsize_12.pdf", height=15, width=22, pointsize=12)
#pdf("15_22_pointsize_18.pdf", height=15, width=22, pointsize=18)
#pdf("15_22_pointsize_24.pdf", height=15, width=22, pointsize=24)
#pdf("15_22_pointsize_30.pdf", height=15, width=22, pointsize=30)

#write log
sessionInfo()
