################################################################################
# Align strep_clades CMC data against 
# S. sinensis genome with the highest mapping from output of 
# bwamap_sinensis_sanguinis_CMC.Snakefile
#
# Irina Velsko 02/02/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/sinensis_sanguinis_mapping/non_ind_plaque/sinensis"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_nonind_plaque/*.gz"):
	SAMPLES[os.path.basename(sample).split(".S")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}-sinensis.mapped_rmdup.bam", sample=SAMPLES.keys()),
        expand("{sample}-sinensis.mapped_rmdup.bam.bai", sample=SAMPLES.keys()),
        expand("{sample}-sinensis.mapped_rmdup.reads.txt", sample=SAMPLES.keys()),
        expand("{sample}-sinensis.mapped_rmdup.cov", sample=SAMPLES.keys()),
        expand("{sample}-sinensis.mapped_rmdup.d.cov", sample=SAMPLES.keys())

rule bwa_aln:
    output:
        temp("{sample}-sinensis.sai")
    message: "Align sample {wildcards.sample} against S. sinensis genome using BWA aln"
    params: 
        reffa = "/mnt/archgen/microbiome_calculus/strep_clades/01-data/Streptococcus_sanguinis_sinensis/ncbi_dataset/data/GCF_000767835.1/GCF_000767835.1_ASM76783v1_genomic.fna",
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
        "{sample}-sinensis.sai"
    output:
        bam = "{sample}-sinensis.mapped.bam",
        reads = temp("{sample}-sinensis.mapped.reads.txt")
    message: "Generate alignment file for sample {wildcards.sample}"
    params:
        reffa = "/mnt/archgen/microbiome_calculus/strep_clades/01-data/Streptococcus_sanguinis_sinensis/ncbi_dataset/data/GCF_000767835.1/GCF_000767835.1_ASM76783v1_genomic.fna",
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
        "{sample}-sinensis.mapped.bam"
    output:
        bam = "{sample}-sinensis.mapped_rmdup.bam",
        reads = "{sample}-sinensis.mapped_rmdup.reads.txt"
    message: "Remove duplicate mapped reads for sample {wildcards.sample}"
    params:
    shell:
        """
        samtools rmdup -s {input} {output.bam}
        bedtools bamtobed -i {output.bam} | cut -f 4 | sort > {output.reads}
        """

rule samtools_index_mem_rmdup:
    input:
        "{sample}-sinensis.mapped_rmdup.bam"
    output:
        "{sample}-sinensis.mapped_rmdup.bam.bai"
    message: "Index bam file for sample {wildcards.sample} mapped against cat'd genes"
    params:
    shell:
        """
        samtools index {input}
        """

rule bedtools_coverage_rmdup:
    input:
        "{sample}-sinensis.mapped_rmdup.bam"
    output:
        "{sample}-sinensis.mapped_rmdup.cov"
    message: "Calculate depth/breadth of coverage of {wildcards.sample}"
    params:
    shell:
        """
        bedtools genomecov -ibam {input} > {output}
        """

rule bed_cov_rmdup_d:
    input:
        "{sample}-sinensis.mapped_rmdup.bam"
    output:
        "{sample}-sinensis.mapped_rmdup.d.cov"
    message: "Calculate depth/breadth of coverage of {wildcards.sample}"
    params:
    shell:
        """
        bedtools genomecov -ibam {input} -d > {output}
        """

