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

###########
# GLOBALS #
###########

fai_file <- snakemake@input[["fai_file"]]
st_depth_file <- snakemake@input[["st_depth_file"]]

########
# MAIN #
########

st_depth_names <- c("Chr", "BP", "depth")
st_depth <- fread(st_depth_file, col.names = st_depth_names)

fai_names <- c("Chr", "chr_length")
##read in fai file, first 2 columns (chr name and length)
fai <- fread(fai_file, select = 1:2, col.names = fai_names)

n_to_label <- 20

# chromosome coordinates
chr_coords <- copy(fai)
##order by chr length (longest to shortest)
setorder(chr_coords, -chr_length, Chr)
##calculate end of chr = length of chromosome + length of all chromosomes beforehand
chr_coords[, chr_end := cumsum(chr_length)]
##calc start = end position - length of chr +1
chr_coords[, chr_start := chr_end - chr_length + 1]
##label position = middle of each chromosome?
chr_coords[, lab_pos := chr_start + round(mean(c(chr_start, chr_end)), 0), by = Chr]
###???
pos_to_label = chr_coords[, seq(1, max(chr_end), length.out = n_to_label)]

label_positions <- sapply(pos_to_label, function(x)
  chr_coords[, .I[which.min(abs(x - lab_pos))]])

chr_coords[label_positions, x_lab := Chr]
chr_coords[is.na(x_lab), x_lab := ""]

# homemade manhattan plot
st_depth_with_len <- merge(st_depth, chr_coords, by="Chr")
setorder(st_depth_with_len, -chr_length, Chr, BP)
st_depth_with_len[, bp_coord := BP + chr_start - 1]

# pick out the outliers
q99 <- st_depth_with_len[, quantile(`depth`, 0.99)]
st_depth_with_len[`depth` > q99, outlier := TRUE]
st_depth_with_len[outlier == TRUE, point_colour := Chr]
st_depth_with_len[is.na(outlier), point_colour := NA]

##skip##
# order the contigs - mixedsort orders strings containing characters and numbers
#st_depth_with_len[, point_colour := factor(point_colour, levels = unique(gtools::mixedsort(point_colour, na.last = TRUE)))]

# include labels for scaffold types

##list of viral scaffolds
full_viral_scaffold_list <- fread("data/viral_scaffold_ids.txt", header=FALSE)
dedup_viral_list <- unique(full_viral_scaffold_list)
##80 and 8624 have non-viral hits on them & were scaffolded, 3939 and 28315 higher depth and non-viral hits
prob_not_viral <- list("Scaffold80", "Scaffold8624", "Scaffold3939", "Scaffold28315")
viral_scaffold_table <- subset(dedup_viral_list, !(V1 %in% prob_not_viral))
##scaffold ids edited for hic mapping
viral_scaffold_table$V1 <- paste(viral_scaffold_table$V1, "unscaffolded", sep="__")
viral_scaffold_list <- list(viral_scaffold_table$V1)
fwrite(viral_scaffold_list, snakemake@output[["viral_scaffold_list"]])

##list of hi-c scaffolds
PGA_depth <- dplyr::filter(st_depth, grepl("PGA_scaffold", Chr))
hic_scaffold_list <- list(unique(PGA_depth$Chr))
fwrite(hic_scaffold_list, snakemake@output[["hic_scaffold_list"]])

##make labels and subset data
st_depth_boxpl <- st_depth_with_len
st_depth_boxpl$label <- ifelse(st_depth_boxpl$Chr %in% viral_scaffold_list, 'viral, unscaffolded',
                               ifelse(st_depth_boxpl$Chr %in% hic_scaffold_list, 'hi-c scaffold', 'unscaffolded'))
st_depth_boxpl <- subset(st_depth_boxpl, !(label == "unscaffolded"))

##################
##manhattan plot##
##################

#pdf(snakemake@output[["manhattan_plot"]], height=30, width=20)
#ggplot() +
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 30,
#                                   hjust = 1,
#                                   vjust = 1),
#        axis.ticks.x = element_blank(),
#        axis.ticks.length.x = unit(0, "mm"),
 #       panel.grid.major.x = element_blank(),
#        panel.grid.minor.x = element_blank(),
#        legend.position = "none") +
#  scale_x_continuous(breaks = chr_coords[, lab_pos],
#                     labels = chr_coords[, x_lab]) +
#  scale_colour_viridis_d() +
#  geom_hline(yintercept = q99) +
#  geom_point(mapping = aes(x = bp_coord,
#                           y = `depth`),
#             data = st_depth_with_len[is.na(point_colour)]) +
#  geom_point(mapping = aes(x = bp_coord,
#                           y = `depth`,
#                           colour = point_colour),
#             data = st_depth_with_len[!is.na(point_colour)])
#dev.off()

##or try with ggsave 
#ggsave("gg_manhattan.pdf", manhattan_plot)

#d99 <- st_depth_with_len[, quantile(`Smoothed D_est`, 0.99)]
#ggplot(st_depth_with_len, aes(x = bp_coord, y = `Smoothed D_est`)) +
# theme(axis.text.x = element_text(angle = 30,
#                                  hjust = 1,
#                                  vjust = 1),
#       axis.ticks.x = element_blank(),
#       axis.ticks.length.x = unit(0, "mm")) +
# scale_x_continuous(breaks = chr_coords[, lab_pos],
#                    labels = chr_coords[, x_lab]) +
# geom_hline(yintercept = d99) +
# geom_point()


#write log
sessionInfo()