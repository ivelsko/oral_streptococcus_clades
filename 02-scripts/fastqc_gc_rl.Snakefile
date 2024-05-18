################################################################################
# Run FastQC on eager-processed reads
#
# Irina Velsko, 20/09/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/sample_gc_rl/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/projects1/microbiome_calculus/iberian/03-preprocessing/screening/eager2/samtools/filter/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################
# /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/sourceTracker_sources/
## /mnt/archgen/microbiome_calculus/Cameroon_plaque/03-preprocessing/all_data_combined/
## /mnt/archgen/microbiome_calculus/Cameroon_plaque/03-preprocessing/all_data_combined/hmp_pq_eager2/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/asangba2022/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/breally2020/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/brito2016/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/eerkens2018/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/fagernaes2021/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/granehaell2021/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/hmp2012_suprapq/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/jacobson2020/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/moraitou2022/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/neukamm2020/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/ottoni2019/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/ozga2019/
## /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/weyrich2017/
## /mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/mann2018_world/eager2_out/
## /mnt/archgen/microbiome_calculus/RIII_simple/03-preprocessing/by_sample/eager2_out/
# /mnt/archgen/projects1/microbiome_calculus/iberian/03-preprocessing/screening/

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input:
        expand("FastQC/{sample}.unmapped_fastqc.zip", sample=SAMPLES.keys())

rule fastqc:
    output:
        "FastQC/{sample}.unmapped_fastqc.zip"
    message: "Run FASTQC on {wildcards.sample}"
    params:
        fasta = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        fastqc {params.fasta} -o /mnt/archgen/microbiome_calculus/strep_clades/03-preprocessing/sample_gc_rl/FastQC
        """
