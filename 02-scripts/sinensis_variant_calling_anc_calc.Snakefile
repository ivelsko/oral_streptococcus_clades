################################################################################
# Call variants on the ancient calculus samples mapped to S. sinensis
#
# Irina Velsko, 02/02/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/sinensis_sanguinis_mapping/ancient_calculus/sinensis"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/sinensis_sanguinis_mapping/ancient_calculus/sinensis/*.mapped_rmdup.bam"):
	SAMPLES[os.path.basename(sample).split("-s")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.sinensis.vcf", sample=SAMPLES.keys()),
        expand("{sample}.sinensis.vcf.tsv", sample=SAMPLES.keys())

rule call_variants:
    output:
        "{sample}.sinensis.vcf",
    message: "Run {wildcards.sample} mapped to S. sinensis through bcftools to call variants"
    params: 
        reffa = "/mnt/archgen/microbiome_calculus/strep_clades/01-data/Streptococcus_sanguinis_sinensis/ncbi_dataset/data/GCF_000767835.1/GCF_000767835.1_ASM76783v1_genomic.fna",
        bam = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        bcftools mpileup -Ou -f {params.reffa} {params.bam} | \
            bcftools call -Ou -mv -A --ploidy 1 | \
            bcftools filter -s LowQual -e '%QUAL<20' > {output}
        """

rule convert_file:
    input:
        "{sample}.sinensis.vcf"
    output:
        "{sample}.sinensis.vcf.tsv"
    message: "Run {wildcards.sample} S. sinensis map vcf to tsv"
    params: 
    shell:
        """
        bcftools norm -m - {input} | \
            bcftools query -f '%CHROM\t%POS\t%REF\t%DP\t%ALT\n' > {output}
        """
