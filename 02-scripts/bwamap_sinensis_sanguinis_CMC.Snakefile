################################################################################
# Align strep_clades ancient calculus data against 
# representative S. sinensis and S. sanguinis genomes
#
# Irina Velsko 10/01/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/sinensis_sanguinis_mapping/non_ind_plaque"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_nonind_plaque/*.gz"):
	SAMPLES[os.path.basename(sample).split(".S")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}-sinensis_sanguinis_rep.mapped_rmdup.bam", sample=SAMPLES.keys()),
        expand("{sample}-sinensis_sanguinis_rep.mapped_rmdup.bam.bai", sample=SAMPLES.keys()),
        expand("{sample}-sinensis_sanguinis_rep.mapped_rmdup.reads.txt", sample=SAMPLES.keys()),
        expand("{sample}-sinensis_sanguinis_rep.mapped_rmdup.cov", sample=SAMPLES.keys())

rule bwa_aln:
    output:
        temp("{sample}-sinensis_sanguinis_rep.sai")
    message: "Align sample {wildcards.sample} against S. sinensis/S. sanguinis references using BWA mem"
    group: "bwa"
    params: 
        reffa = "/mnt/archgen/microbiome_calculus/strep_clades/01-data/Streptococcus_sanguinis_sinensis/rep_genomes/S_sinensis_sanguinis.fna",
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 4
    shell:
        """
        bwa aln -l 32 -n 0.01 -t {threads} \
            {params.reffa} \
            {params.fastq} > {output}
        """

rule bwa_samse_aln:
    input:
        "{sample}-sinensis_sanguinis_rep.sai"
    output:
        bam = "{sample}-sinensis_sanguinis_rep.mapped.bam",
        reads = temp("{sample}-sinensis_sanguinis_rep.mapped.reads.txt")
    message: "Generate alignment file for sample {wildcards.sample}"
    group: "bwa"
    params:
        reffa = "/mnt/archgen/microbiome_calculus/strep_clades/01-data/Streptococcus_sanguinis_sinensis/rep_genomes/S_sinensis_sanguinis.fna",
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        bwa samse \
            {params.reffa} \
            {input} \
            {params.fastq} | \
        samtools view -Sb -F 4 - | samtools sort - -o {output.bam} 
        bedtools bamtobed -i {output.bam} | cut -f 4 | sort > {output.reads}
        """

rule samtools_rmdup_mem:
    input:
        "{sample}-sinensis_sanguinis_rep.mapped.bam"
    output:
        bam = "{sample}-sinensis_sanguinis_rep.mapped_rmdup.bam",
        reads = "{sample}-sinensis_sanguinis_rep.mapped_rmdup.reads.txt"
    message: "Remove duplicate mapped reads for sample {wildcards.sample}"
    group: "bwa"
    params:
    shell:
        """
        samtools rmdup -s {input} {output.bam}
        bedtools bamtobed -i {output.bam} | cut -f 4 | sort > {output.reads}
        """

rule samtools_index_mem_rmdup:
    input:
        "{sample}-sinensis_sanguinis_rep.mapped_rmdup.bam"
    output:
        "{sample}-sinensis_sanguinis_rep.mapped_rmdup.bam.bai"
    message: "Index bam file for sample {wildcards.sample} mapped against cat'd genes"
    group: "bwa"
    params:
    shell:
        """
        samtools index {input}
        """

rule bedtools_coverage_rmdup:
    input:
        "{sample}-sinensis_sanguinis_rep.mapped_rmdup.bam"
    output:
        "{sample}-sinensis_sanguinis_rep.mapped_rmdup.cov"
    message: "Calculate depth/breadth of coverage of {wildcards.sample}"
    group: "bwa"
    params:
    shell:
        """
        bedtools genomecov -ibam {input} > {output}
        """

