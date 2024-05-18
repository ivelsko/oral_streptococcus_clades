################################################################################
# Combine the MPA-style report files from Kraken2 wGTDB r202
#
# Irina Velsko, 01/06/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/output_mod_calc/*.report_mpa.tsv"):
	SAMPLES[os.path.basename(sample).split("_")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        "output/{sample}.kraken2.gtdb.report_mpa.tsv"

rule combine_mpa:
    output:
        "combined.mod_calc.kraken2.gtdb.report_mpa.tsv",
    message: "Run {wildcards.sample} through kraken2 database w/bacteria, archaea, human"
    params: 
        infiles = expand("{sample}", sample=SAMPLES.values()),
    threads: 32
    shell:
        """
        /mnt/archgen/users/velsko/bin/krakenTools/combine_mpa.py \
        -i {params.infiles}  \
        -o {output}
        """

# rule combine_mpa:
#     output:
#         "combined.mod_calc.kraken2.gtdb.report_mpa.tsv",
#     message: "Run {wildcards.sample} through kraken2 database w/bacteria, archaea, human"
#     params: 
#         infiles = expand("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/output_mod_calc/*.report_mpa.tsv")
#     threads: 32
#     shell:
#         """
#         /mnt/archgen/users/velsko/bin/krakenTools/combine_mpa.py \
#         -i {params.infiles}  \
#         -o {output}
#         """
