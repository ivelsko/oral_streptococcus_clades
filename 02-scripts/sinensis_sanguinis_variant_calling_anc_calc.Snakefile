################################################################################
# Call variants on the ancient calculus samples mapped to S. sinensis & S. sanguinis
#
# Irina Velsko, 17/01/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/sinensis_sanguinis_mapping/ancient_calculus/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/sinensis_sanguinis_mapping/ancient_calculus/*.mapped_rmdup.bam"):
	SAMPLES[os.path.basename(sample).split("-s")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.sinensis_sanguinis.vcf", sample=SAMPLES.keys())

rule call_variants:
    output:
        "{sample}.sinensis_sanguinis.vcf",
    message: "Run {wildcards.sample} mapped to S. sinensis and S. sanguinis through bcftools to call variants"
    params: 
        reffa = "/mnt/archgen/microbiome_calculus/strep_clades/01-data/Streptococcus_sanguinis_sinensis/rep_genomes/S_sinensis_sanguinis.fna",
        bam = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        bcftools mpileup -Ou -f {params.reffa} {params.bam} | \
            bcftools call -Ou -mv --ploidy 1 | \
            bcftools filter -s LowQual -e 'QUAL<20' > {output}
        """

#            bcftools filter -s LowQual -e '%QUAL<20 || DP>2' > {output}
