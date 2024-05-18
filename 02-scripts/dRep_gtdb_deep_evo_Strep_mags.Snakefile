################################################################################
# Run dRep on Strep genomes of unknown and known Strep groups
# to place the unknowns in groups based on ANI, if possible
#
# Irina Velsko, 19/05/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/deep_evo_mag_drep"


if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        "dRep_out/data_tables/Ndb.csv"
#         "dRep_out/data_tables/Ndb.csv"

rule dRep:
    output:
        "dRep_out/data_tables/Ndb.csv"
    message: "Run dRep to cluster Strep genomes for gtdb unknown species"
    params: 
        genomes = "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/deep_evo_mag_drep/strep_genomes_list.tsv"
    threads: 24
    shell:
        """
        dRep compare -p {threads} dRep_out/ -g {params.genomes} --S_algorithm ANImf -pa 0.95 -sa 0.99
        """
     
        