################################################################################
# Call variants on the CMC samples mapped to S. sanguinis
#
# Irina Velsko, 02/02/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/sinensis_sanguinis_mapping/non_ind_plaque/sanguinis"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/sinensis_sanguinis_mapping/non_ind_plaque/sanguinis/*.mapped_rmdup.bam"):
	SAMPLES[os.path.basename(sample).split("-")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.sanguinis.vcf", sample=SAMPLES.keys()),
        expand("{sample}.sanguinis.vcf.tsv", sample=SAMPLES.keys())

rule call_variants:
    output:
        "{sample}.sanguinis.vcf",
    message: "Run {wildcards.sample} mapped to S. sanguinis through bcftools to call variants"
    params: 
        reffa = "/mnt/archgen/microbiome_calculus/strep_clades/01-data/Streptococcus_sanguinis_sinensis/ncbi_dataset/data/GCF_013343115.1/GCF_013343115.1_ASM1334311v1_genomic.fna",
        bam = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        bcftools mpileup -Ou -f {params.reffa} {params.bam} | \
            bcftools call -Ou -mv -A --ploidy 1 | \
            bcftools filter -s LowQual -e 'QUAL<20' > {output}
        """

rule convert_file:
    input:
        "{sample}.sanguinis.vcf"
    output:
        "{sample}.sanguinis.vcf.tsv"
    message: "Run {wildcards.sample} S. sanguinis map vcf to tsv"
    params: 
    shell:
        """
        bcftools norm -m - {input} | \
            bcftools query -f '%CHROM\t%POS\t%REF\t%DP\t%ALT\n' > {output}
        """
