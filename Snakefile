#!/usr/bin/env python3

import pathlib2
import pandas
import os

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

####for bam files####

def find_bam_files(bam_dir):
#Make list of files
    path_generator = os.walk(peptide_dir, followlinks = True)
    my_files = list((dirpath, filenames)
        for (dirpath, dirname, filenames)
        in path_generator)
#Make new dictionary & populate with files
    my_peptide_files = {}
    for dirpath, filenames in my_files:
        for filename in filenames:
            if filename.endswith('.bam'):
                my_peptide_files = str(pathlib2.Path(dirpath))
    return(my_peptide_files)

###########
# GLOBALS #
###########

####for bam files####
bam_sample_key_file = 'data/bam_sample_key.csv'
bam_sample_key = pandas.read_csv(bam_sample_key_file)
bam_dir = 'data/hyp_bams'
all_samples = sorted(set(bam_sample_key['filename']))
contig_id_key_file = 'data/viral_hic_contig_ids.csv'
contig_id_key = pandas.read_csv(contig_id_key_file)
all_contigs = sorted(set(contig_id_key))

##############
# CONTAINERS #
##############

bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'

#########
# SETUP #
#########

# generate name to filename dictionary
all_samples = sorted(set(bam_sample_key['filename']))

#########
# RULES #
#########

rule target:
    input:
        expand('output/mean_depth/{sample}_mean_depth_table.csv', sample=all_samples),
        'data/Mh_Hi-C_PGA_assembly.fasta.fai',
        'output/viral_scaffold_list.txt',
        expand('output/depth_analysis/{sample}_boxplot.pdf', sample=all_samples),
        'output/depth_analysis/boxplot_panel.pdf',
        #'output/depth_analysis/indiv12_manhattan.pdf',
        expand('output/samtools_depth/filtered/filtered_{sample}_depth.out', sample=all_samples),
        expand('output/depth_analysis/{sample}_viral_hic_depth_wilcox.txt', sample=all_samples),
        expand('output/depth_analysis/{sample}_viral_hist.pdf', sample=all_samples),
        expand('output/depth_analysis/{sample}_hic_hist.pdf', sample=all_samples),
        expand('output/depth_analysis/kruskal_wallis/{sample}_kw_test.txt', sample=all_samples)

rule depth_kruskal_wallis:
    input:
        st_depth_file = 'output/samtools_depth/filtered/filtered_{sample}_depth.out',
        scaffold_id_table = 'data/scaffold_id_table.csv'
    output:
        kw_res = 'output/depth_analysis/kruskal_wallis/{sample}_kw_test.txt',
        pw_res = 'output/depth_analysis/kruskal_wallis/{sample}_pw_test.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/kruskal_wallis/{sample}.log'
    script:
        'src/depth_kruskal_wallis.R'

rule depth_t_test:
    input:
        st_depth_file = 'output/samtools_depth/filtered/filtered_{sample}_depth.out',
        scaffold_id_table = 'data/scaffold_id_table.csv'
    output:
        unpaired_ttest_results = 'output/depth_analysis/{sample}_viral_hic_depth_unpaired_ttest.txt',
        viral_depth_hist = 'output/depth_analysis/{sample}_viral_hist.pdf',
        hic_depth_hist = 'output/depth_analysis/{sample}_hic_hist.pdf',
        wilcox_res = 'output/depth_analysis/{sample}_viral_hic_depth_wilcox.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/ttests/{sample}.log'
    script:
        'src/depth_ttest.R'

rule depth_boxplot_panel:
    input:
        st_depth_file_i12 = 'output/samtools_depth/filtered/filtered_indiv12_mhyp_lincoln_depth.out',
        st_depth_file_i28 = 'output/samtools_depth/filtered/filtered_indiv28_mhyp_lincoln_depth.out',
        st_depth_file_i52 = 'output/samtools_depth/filtered/filtered_indiv52_mhyp_lincoln_depth.out',
        st_depth_file_i84 = 'output/samtools_depth/filtered/filtered_indiv84_mhyp_lincoln_depth.out',
        scaffold_id_table = 'data/scaffold_id_table.csv'
    output:
        boxplot_panel = 'output/depth_analysis/boxplot_panel.pdf'
    singularity:
        tidyverse_container
    log:
        'output/logs/boxplots/panel_depth_boxplot.log'
    script:
        'src/depth_boxplot_panel.R'

rule depth_boxplot:
    input:
        samtools_depth = 'output/samtools_depth/filtered/filtered_{sample}_depth.out',
        scaffold_id_table = 'data/scaffold_id_table.csv'
    output:
        boxplot_y_zoom = 'output/depth_analysis/{sample}_boxplot_y_zoom.pdf',
        boxplot = 'output/depth_analysis/{sample}_boxplot.pdf'
    singularity:
        tidyverse_container
    log:
        'output/logs/boxplots/{sample}_depth_boxplot.log'
    script:
        'src/depth_boxplot.R'

##don't think this works?
rule manhattan_plot:
    input:
        fai_file = 'data/Mh_Hi-C_PGA_assembly.fasta.fai',
        st_depth_file = 'output/samtools_depth/filtered/filtered_indiv12_mhyp_lincoln_depth.out'
    output:
        viral_scaffold_list = 'output/viral_scaffold_list.txt',
        hic_scaffold_list = 'output/hic_scaffod_list.txt',
        #manhattan_plot = 'output/depth_analysis/indiv12_manhattan.pdf'
    singularity:
        tidyverse_container
    threads:
        20
    log:
        'output/logs/manhattan_plot.log'
    script:
        'src/manhattan_plot.R'

rule hic_genome_fai:
    input:
        hic_genome = 'data/Mh_Hi-C_PGA_assembly.fasta'
    output:
        hic_fai = 'data/Mh_Hi-C_PGA_assembly.fasta.fai'
    singularity:
        samtools_container
    threads:
        10
    log:
        'output/logs/hic_genome_fai.log'
    shell:
        'samtools faidx '
        '{input.hic_genome}'

##calc read depth across scaffolds
rule calc_mean_depth:
     input:  
         depth = 'output/samtools_depth/{sample}_depth.out'
     output:
         mean_depth_table = 'output/mean_depth/{sample}_mean_depth_table.csv'
     singularity:
         tidyverse_container
     threads:
         10
     log:
         'output/logs/calc_mean_depth/{sample}.log'
     script:
         'src/calc_mean_scaffold_depth.R'

##can I loop through this separating each sample depth file into chromosome files somehow?
rule filter_depth_file:
    input:
        depth_out = 'output/samtools_depth/{sample}_depth.out',
        viral_hic_ids_list = 'output/viral_and_hic_scaffold_ids.txt'
    output:
        filtered_depth = 'output/samtools_depth/filtered/filtered_{sample}_depth.out'
    shell:
        'egrep -wf {input.viral_hic_ids_list} {input.depth_out} > {output.filtered_depth}'

rule samtools_depth:
    input:
        sorted_bam = 'data/hyp_bams/{sample}.bam'
    output:
        depth_out = 'output/samtools_depth/{sample}_depth.out'
    log:
        'output/logs/samtools_depth/{sample}.log'
    threads:
        20
    singularity:
        samtools_container
    shell:
        'samtools depth '
        '{input.sorted_bam} '
        '-a '
        '> {output.depth_out} '
        '2> {log}'