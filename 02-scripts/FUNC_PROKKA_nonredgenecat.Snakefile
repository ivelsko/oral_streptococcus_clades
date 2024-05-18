####################################################################################################
# Project: Strep clade distribution in oral samples
# Part: Evaluating the functional capacity
# Step: Construction and functional profiling of a non-redundant gene catalogue
#       from de-novo assembled contigs using PROKKA annotations
#
#
# Alex Huebner
# modified by Irina Velsko for Strep clades
####################################################################################################

import bz2
from glob import glob
import os
import re

import pandas as pd
import pyfastx

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_anc_calc/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################


#### SAMPLES ###################################################################
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_anc_calc
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_nhp_baboon
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_nhp_chimp
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_nhp_general_site
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_nhp_gorilla
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_nhp_howler
###
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_mod_calc/*.gz
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_ind_plaque/*.gz
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_nonind_plaque/*.gz
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_ind_buccal_mucosa/*.gz
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_nonind_buccal_mucosa/*.gz
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_ind_saliva/*.gz
# /mnt/archgen/microbiome_calculus/strep_clades/04-analysis/kraken2_gtdb_r202/input_nonind_saliva/*.gz
################################################################################


rule all:
    input:
        expand("99-analysis/tmp/nonredgenecat/{sample}.gmgc10_95nr_complete.sam.bz2", sample=SAMPLES.keys()),
        expand("99-analysis/tmp/nonredgenecat/{sample}.gmgc10_95nr_complete.sorted.bam", sample=SAMPLES.keys())

#### Alignment against GMGC ########################################################################

rule download_gmgc_catalogue:
    output:
        "/mnt/archgen/users/huebner/refdbs/GMGC10.data/subcatalogs/GMGC10.95nr.complete.fna.gz"
    message: "Download the Global Microbial Gene Catalogue (all complete genes; 95% redundancy; only complete genes)"
    params:
        url = "git@git.embl.de:coelho/GMGC10.data.git",
        dir = "99-analysis/gmgc",
        subcat = "subcatalogs/GMGC10.95nr.complete.fna.gz"
    shell:
        """
        cd {params.dir}
        git clone {params.url} 
        git-annex get {params.subcat}
        """

rule build_bowtie2_index:
    input:
        "/mnt/archgen/users/huebner/refdbs/GMGC10.data/subcatalogs/GMGC10.95nr.complete.fna.gz"
    output:
        "/mnt/archgen/users/huebner/refdbs/GMGC10.data/subcatalogs/GMGC10.95nr.complete.1.bt2l"
    message: "Build index of GMGC for BowTie2"
    conda: "MISC_BowTie2.yaml"
    resources:
        mem = 200
    params:
        index = "/mnt/archgen/users/huebner/refdbs/GMGC10.data/subcatalogs//GMGC10.95nr.complete"
    threads: 16
    shell:
        """
        bowtie2-build --threads {threads} -f {input} {params.index}
        """

rule bowtie2_alignment:
    input:
        "/mnt/archgen/users/huebner/refdbs/GMGC10.data/subcatalogs/GMGC10.95nr.complete.1.bt2l"
    output:
        "99-analysis/tmp/nonredgenecat/{sample}.gmgc10_95nr_complete.sam.bz2"
    message: "Align against gene catalogue and sort BAM file: {wildcards.sample}"
    conda: "MISC_BowTie2.yaml"
    resources:
        mem = 180,
        cores = 16
    params:
        db = "/mnt/archgen/users/huebner/refdbs/GMGC10.data/subcatalogs/GMGC10.95nr.complete",
        fastq = lambda wildcards: SAMPLES[wildcards.sample],
        nmismatches = lambda wildcards: 1 if wildcards.sample[1:3] == "SM" else 0
    benchmark: "snakemake_tmp/bowtie2_alignment.{sample}.benchmark"
    threads: 16
    shell:
        """
        bowtie2 -x {params.db} \
                -U {params.fastq} \
                -D 20 -R 3 -N {params.nmismatches} -L 20 -i S,1,0.50 \
                --no-unal \
                --threads {threads} | bzip2 > {output}
        """

rule sam2bam:
    input:
        "99-analysis/tmp/nonredgenecat/{sample}.gmgc10_95nr_complete.sam.bz2"
    output:
        "99-analysis/tmp/nonredgenecat/{sample}.gmgc10_95nr_complete.reheader.sam.bz2"
    message: "Extract genes to which reads were aligned and convert to BAM: {wildcards.sample}"
    run:
        # Due to the too many genes in the gene catalogue, we cannot directly
        # create a BAM file because there is limit of the number of genome
        # entries in a BAM header. Therefore, we will create a list of all
        # genes to which reads were aligned and then adapt the header
        gene_ids = {}
        genes = []
        with bz2.open(input[0], "rt") as samfile:
            gene_id = 0
            for line in samfile:
                if line.startswith("@SQ"):
                    gene_ids[line.split("\t")[1].split(":")[1]] = gene_id
                    gene_id += 1
                elif not line.startswith("@"):
                    gene = line.rstrip().split("\t")[2]
                    genes.append(gene_ids[gene])
        found_genes = set(genes)

        with bz2.open(output[0], "wt") as outfile:
            with bz2.open(input[0], "rt") as samfile:
                sq_id = 0
                for line in samfile:
                    if line.startswith("@SQ"):
                        if sq_id in found_genes:
                            outfile.write(line)
                        sq_id += 1
                    else:
                        outfile.write(line)

rule sort_bam:
    input:
        "99-analysis/tmp/nonredgenecat/{sample}.gmgc10_95nr_complete.reheader.sam.bz2"
    output:
        "99-analysis/tmp/nonredgenecat/{sample}.gmgc10_95nr_complete.sorted.bam"
    conda: "FUNC_PROKKA_nonredgenecat.yaml"
    resources:
        mem = 16,
        cores = 3
    threads: 3
    shell:
        """
        bzcat {input} | \
        samtools view -Sb - | \
        samtools sort -T {wildcards.sample} -o {output} -
        """