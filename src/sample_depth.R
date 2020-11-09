library(data.table)
library(readr)
library(dplyr)
library(ggplot2)

##find paths for all mean depth tables
mean_depth_tables <- list.files("output/mean_depth",
                                pattern="_mean_depth_table.csv", full.names=TRUE)
##read in list of all results
depth_results_list <- lapply(mean_depth_tables, fread, skip=1)
##name results based of filename
names(depth_results_list) <- gsub("output/mean_depth/(.+)_mean_depth_table.csv", "\\1", mean_depth_tables)
##bind rows in list based off filename - makes full table
all_depth <- dplyr::bind_rows(depth_results_list, .id="filename")
setnames(all_depth, old=c("V1", "V2"), new=c("scaffold_id", "mean_depth"))
##generate list of all samples
sample_ids <- unique(all_depth$`filename`)

##take depth table, and calculate mean for all scaffolds belonging to each sample in table
MEAN_DEPTH <- function(x, depth){
  sample_scaffolds<-depth[depth$`filename`==x,]
  mean_depth<-mean(sample_scaffolds$mean_depth)
  return(data.table(sample_id=x, mean_depth=mean_depth))
}

##################
##Hi-C Scaffolds##
##################

##PGA scaffold depth
PGA_depth <- dplyr::filter(all_depth, grepl("PGA_scaffold", scaffold_id))
sample_PGA_depths <- lapply(sample_ids, MEAN_DEPTH, depth=PGA_depth)
sample_PGA_depths_table <- rbindlist(sample_PGA_depths)

##unanchored scaffold depth
unscaffolded_depth <- dplyr::filter(all_depth, grepl("unscaffolded", scaffold_id))
sample_unscaffolded_depths <- lapply(sample_ids, MEAN_DEPTH, depth=unscaffolded_depth)
sample_unscaffolded_depths_table <- rbindlist(sample_unscaffolded_depths)

###################
##viral scaffolds##
###################

##need list of viral scaffolds
##need to get better list than this really
##read in list of scaffolds
full_viral_scaffold_list <- fread("data/viral_scaffold_ids.txt", header=FALSE)
dedup_viral_list <- unique(full_viral_scaffold_list)
##80 and 8624 have non-viral hits on them
prob_not_viral <- list("Scaffold80", "Scaffold8624", "Scaffold3939", "Scaffold28315")
viral_scaffold_list <- subset(dedup_viral_list, !(V1 %in% prob_not_viral))
##scaffold ids edited for hic mapping
viral_scaffold_list$V1 <- paste(viral_scaffold_list$V1, "unscaffolded", sep="__")

##subset depth table to get results for viral scaffolds only
viral_depth_table <- subset(all_depth, (all_depth$scaffold_id %in% viral_scaffold_list$V1))
##mean depth of all viral scaffolds per sample
sample_viral_depths <- lapply(sample_ids, MEAN_DEPTH, depth=viral_depth_table)
sample_viral_depths_table <- rbindlist(sample_viral_depths)


##table of PGA mean depth and viral mean depth for each sample
PGA_vs_viral_depth <- merge(sample_PGA_depths_table, sample_viral_depths_table, by="sample_id")
PGA_vs_viral_depth <- merge(PGA_vs_viral_depth ,sample_unscaffolded_depths_table, by="sample_id")
setnames(PGA_vs_viral_depth, old=c("mean_depth.x", "mean_depth.y", "mean_depth"), new=c("PGA_mean_depth", "viral_mean_depth", "unscaff_mean_depth"))
##make row of PGA depth/viral depth for each sample
PGA_vs_viral_depth$PGA_vs_viral <- ((PGA_vs_viral_depth$PGA_mean_depth)/(PGA_vs_viral_depth$viral_mean_depth))
##add column for wasp location
PGA_vs_viral_depth$location <- PGA_vs_viral_depth$sample_id
PGA_vs_viral_depth$location <- tstrsplit(PGA_vs_viral_depth$location, "mhyp_", keep=c(2))
fwrite(PGA_vs_viral_depth, "output/depth_analysis/scaff_unscaff_viral_depth.csv")

##merge with stage info
sample_key <- fread("data/bam_sample_key.csv")
ok_sample_key <- subset(sample_key, (sample_key$Issue=="N"))
plot_data <- merge(PGA_vs_viral_depth, ok_sample_key, by.x="sample_id", by.y="filename")

viral_scaffold_list <- subset(dedup_viral_list, !(V1 %in% prob_not_viral))


##dot plot
ggplot(plot_data, aes(x=PGA_mean_depth, y=viral_mean_depth, colour=Stage, shape=factor(Location)))+
  geom_abline(slope=1, intercept=0,  color="grey48", linetype="dashed")+
  geom_point()+labs(x="Mean Depth for Hi-C Scaffolds", y="Mean Depth for Viral Contigs")+
  coord_fixed(xlim=c(0, 25), ylim=c(0, 25))

##box plot - may need to manipulate table more to plot this with box plot
ggplot(PGA_vs_viral_depth, aes(x=PGA_mean_depth, y=viral_mean_depth, colour=location))+
  geom_boxplot()+labs(x="Mean Depth for Hi-C Scaffolds", y="Mean Depth for Viral Contigs")+
  coord_fixed(xlim=c(0, 25), ylim=c(0, 25))


##compare 2 genome contigs
contig1_depth <- subset(all_depth, (all_depth$scaffold_id=="PGA_scaffold0__1249_contigs__length_14794116"))
contig2_depth <- subset(all_depth, (all_depth$scaffold_id=="PGA_scaffold1__1124_contigs__length_11802476"))
plot_depth_table <- merge(contig1_depth, contig2_depth, by="filename")
##add in location
plot_depth_table_info <- merge(plot_depth_table, ok_sample_key, by="filename")
##plot depth hic chr0 vs chr1
ggplot(plot_depth_table_info, aes(x=mean_depth.x, y=mean_depth.y, colour=Stage, shape=factor(Location)))+
  geom_abline(slope=1, intercept=0,  color="grey48", linetype="dashed")+
  geom_point()+labs(x="Mean Depth for Hi-C Chr 0", y="Mean Depth for Hi-C Chr 1")+
  coord_fixed(xlim=c(0, 25), ylim=c(0, 25))

##read in fai to get table of chromosomes vs length
fai_names <- c("Chr", "chr_length")
##read in fai file, first 2 columns (chr name and length)
fai <- fread("data/Mh_Hi-C_PGA_assembly.fasta.fai", select = 1:2, col.names = fai_names)
##subset into hic and viral tables
hic_fai <- dplyr::filter(fai, grepl("PGA_scaffold", Chr))
viral_fai <- subset(fai, (fai$Chr %in% viral_scaffold_list$V1))


##mean scaffold depth
mean(PGA_vs_viral_depth$PGA_mean_depth)
##mean viral contig depth
mean(PGA_vs_viral_depth$viral_mean_depth)
##of all unscaff
mean(PGA_vs_viral_depth$unscaff_mean_depth)

scaffold_ids <- data.table(unique(all_depth$scaffold_id))
##3683 total scaffolds
##12 = Hi-C scaffolds
##19 = viral contigs
##3652 other unanchored contigs

