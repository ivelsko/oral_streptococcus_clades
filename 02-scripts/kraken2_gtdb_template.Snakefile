################################################################################
# Run Kraken2 with GTDB r202 database on XXXXXXXX samples
#
# Irina Velsko, 29/05/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2/input/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("output/{sample}.kraken2.gtdb.output.tsv", sample=SAMPLES.keys()),
        expand("output/{sample}.kraken2.gtdb.report_mpa.tsv", sample=SAMPLES.keys())

rule kraken_classify:
    output:
        outfmt = "output/{sample}.kraken2.gtdb.output.tsv",
        repfmt = "output/{sample}.kraken2.gtdb.report_mpa.tsv"
    message: "Run {wildcards.sample} through kraken2 database w/bacteria, archaea, human"
    params: 
        database = "/mnt/archgen/microbiome_sciences/reference_databases/built/struo2_GTDB_release202/kraken2",
        infile = lambda wildcards: SAMPLES[wildcards.sample],
    threads: 32
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate kraken2
        set -u

        kraken2 --db {params.database} \
        {params.infile} \
        --threads {threads} \
        --output {output.outfmt} \
        --report {output.repfmt} \
        --report-zero-counts \
        --use-mpa-style \
        --use-names  \
        --gzip-compressed
        """

