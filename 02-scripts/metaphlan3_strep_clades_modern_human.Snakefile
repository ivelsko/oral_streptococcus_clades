################################################################################
# Run MetaPhlAn3 on Pacific calculus samples
#
# Irina Velsko, 07/06/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/metaphlan3"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/metaphlan3/input_modern_human/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb/input_modern_human/*.gz

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.metaphlan3.tsv", sample=SAMPLES.keys())

rule metaphlan3:
    output:
        "{sample}.metaphlan3.tsv"
    message: "Run metaphlan3 {wildcards.sample}"
    params: 
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 12
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate mpa3
        set -u

        metaphlan {params.fastq} --input_type fastq --nproc 12 > {output}
        """
