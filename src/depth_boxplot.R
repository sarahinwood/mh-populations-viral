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

samtools_depth <- snakemake@input[["samtools_depth"]]
##scaffold labels
scaffold_id_table <- snakemake@input[["scaffold_id_table"]]

########
# MAIN #
########

##scaffold ID table
scaffold_table <- fread(scaffold_id_table, header=TRUE)
scaffold_table$scaffold_number <- tstrsplit(scaffold_table$'#Name', "_", keep=c(2))
scaffold_table$scaffold_number <- as.numeric(as.character(scaffold_table$scaffold_number))
setorder(scaffold_table, scaffold_number)
scaffold_table$plot_group <- tstrsplit(scaffold_table$plot_label, " and", keep=c(1))


##depth table
st_depth_names <- c("#Name", "BP", "depth")
st_depth <- fread(samtools_depth, col.names=st_depth_names)

##merge with table of scaffold ids for plotting
st_depth_labels <- merge(st_depth, scaffold_table, by='#Name', all.x=TRUE)
st_depth_boxpl <- subset(st_depth_labels, !(plot_group=="Other contig"))
mh_depth_outliers <- list("scaffold_90", "scaffold_995")
st_depth_boxpl <- subset(st_depth_labels, !(`#Name` %in% mh_depth_outliers))

##order legend labels
st_depth_boxpl$plot_group <- factor(st_depth_boxpl$plot_group, levels=c("Hi-C scaffold", "Viral contig"))

#####################
## making box plot ##
#####################

pdf(snakemake@output[["boxplot_y_zoom"]])
ggplot(st_depth_boxpl, aes(x=reorder(`#Name`, scaffold_number), y=depth, colour=plot_group))+
  geom_boxplot(outlier.shape=NA)+
  theme_bw(base_size=18)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank())+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun=mean, geom="point", colour="grey35")+
  scale_colour_viridis(discrete=TRUE)+
  coord_cartesian(ylim = c(0,35))
dev.off()

jpeg(snakemake@output[["boxplot"]])
ggplot(st_depth_boxpl, aes(x=reorder(`#Name`, scaffold_number), y=depth, colour=plot_group))+
  geom_boxplot()+
  theme_bw(base_size=18)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank())+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun=mean, geom="point", colour="grey35")+
  scale_colour_viridis(discrete=TRUE)
dev.off()

#write log
sessionInfo()
