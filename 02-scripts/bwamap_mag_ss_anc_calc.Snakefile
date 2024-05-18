################################################################################
# Align strep_clades ancient calculus data against DeepEvo Strep mags in the #
# Sanguinis clade and representative S. sinensis and S. sanguinis genomes
#
# Irina Velsko 01/16/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/sinensis_sanguinis_mapping/deep_evo_mags/ancient_calculus"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_anc_calc/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}-mag_ss.mapped.bam", sample=SAMPLES.keys()),
        expand("{sample}-mag_ss.mapped.bam.bai", sample=SAMPLES.keys()),
        expand("{sample}-mag_ss.mapped_rmdup.bam", sample=SAMPLES.keys()),
        expand("{sample}-mag_ss.mapped_rmdup.bam.bai", sample=SAMPLES.keys()),
        expand("{sample}-mag_ss.mapped_rmdup.reads.txt", sample=SAMPLES.keys()),
        expand("{sample}-mag_ss.mapped_rmdup.cov", sample=SAMPLES.keys())

rule bwa_aln:
    output:
        temp("{sample}-mag_ss.sai")
    message: "Align sample {wildcards.sample} against S. sinensis/S. sanguinis references using BWA aln"
    params: 
        reffa = "/mnt/archgen/microbiome_calculus/strep_clades/01-data/Streptococcus_sanguinis_sinensis/with_deep_evo_mags/mag_ss_strep.fa",
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
        "{sample}-mag_ss.sai"
    output:
        bam = "{sample}-mag_ss.mapped.bam",
        reads = temp("{sample}-mag_ss.mapped.reads.txt")
    message: "Generate alignment file for sample {wildcards.sample}"
    params:
        reffa = "/mnt/archgen/microbiome_calculus/strep_clades/01-data/Streptococcus_sanguinis_sinensis/with_deep_evo_mags/mag_ss_strep.fa",
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
        "{sample}-mag_ss.mapped.bam"
    output:
        bam = "{sample}-mag_ss.mapped_rmdup.bam",
        reads = "{sample}-mag_ss.mapped_rmdup.reads.txt"
    message: "Remove duplicate mapped reads for sample {wildcards.sample} against cat'd genomes"
    params:
    shell:
        """
        samtools rmdup -s {input} {output.bam}
        bedtools bamtobed -i {output.bam} | cut -f 4 | sort > {output.reads}
        """

rule samtools_index_mem:
    input:
        "{sample}-mag_ss.mapped.bam"
    output:
        "{sample}-mag_ss.mapped.bam.bai"
    message: "Index bam file for sample {wildcards.sample} mapped against cat'd genomes"
    params:
    shell:
        """
        samtools index {input}
        """

rule samtools_index_mem_rmdup:
    input:
        "{sample}-mag_ss.mapped_rmdup.bam"
    output:
        "{sample}-mag_ss.mapped_rmdup.bam.bai"
    message: "Index bam file for sample {wildcards.sample} mapped against cat'd genomes"
    params:
    shell:
        """
        samtools index {input}
        """

rule bedtools_coverage_rmdup:
    input:
        "{sample}-mag_ss.mapped_rmdup.bam"
    output:
        "{sample}-mag_ss.mapped_rmdup.cov"
    message: "Calculate depth/breadth of coverage of {wildcards.sample} mapped against cat'd genomes"
    params:
    shell:
        """
        bedtools genomecov -ibam {input} > {output}
        """

