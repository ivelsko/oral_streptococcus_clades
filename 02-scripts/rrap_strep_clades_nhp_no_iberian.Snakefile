################################################################################
# Run RRAP on strep_clade non-human primate samples w/o Iberian dataset
#
# Irina Velsko, 10/02/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/gmgc_mapping/rrap"

#### SAMPLES ###################################################################
# SAMPLES = {}
# for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/metaphlan3/input_nhp/*.gz"):
# 	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        "rrap_nhp_no_iberian.done"

rule rrap:
    output:
        touch("rrap_nhp_no_iberian.done")
    message: "Run run RRAP on non-human primate samples"
    params: 
        fastq = "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/gmgc_mapping/nhp_input_path_no_iberian.tsv",
        catref = "/mnt/archgen/users/huebner/refdbs/GMGC10.data/subcatalogs/GMGC10.95nr.complete.fna.gz"
    threads: 12
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate rrap_env
        set -u

        rrap -i {params.fastq} \
             -rg {params.catref} \
             -o /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/gmgc_mapping/rrap/output_nhp_no_iberian \
             -n GMGC10.95nr.complete \
             -suffix _1.fastq.gz \
             --skip-indexing \
             --skip-rr \
             --threads {threads}\
             -q
        """
