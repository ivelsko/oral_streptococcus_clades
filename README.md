# Oral Streptococcus Clades
This repo contains all the data tables and scripts needed to repeat the analyses for the manuscript "TBA".

This repository is for keeping the data processing and analysis of oral
microbiome samples from globally and temporally diverse human and non-human
primate populations. This includes taxonomic and functional analyses. The
abundance and distribution of _Streptococcus_ species within and between groups is
assessed. 

Overview about the repository  
`00-documentation`: files with information about samples or _Streptococcus_ species that are used as input for scripts in `02-scripts`.  
`02-scripts`: scripts for sample processing and analysis. Rmd files contain data analysis and supplemental figure scripts.  
`05-results`: results from the data analysis that can stand alone, meaning they aren't strictly used as input for another step.  
`06-publication`: publicaiton-asscociated files for figures and tables  

Metadata for all samples included in the study are in `00-documentation/strep_clades_metadata.tsv`. Submission scripts and notes for all analyses are in `00-documentation/strep_clades_notes.tsv`.

Taxonomic tables for all samples are in `05-results`. There are separate tables
for each sample type. The tables for ancient dental calculus and chipmpanze
dental calculus have sample names for column headers, all other tables have
"Sample #1", etc as headers. To get the sample names, use the corresponding
`<sample_type>_names.tsv` files. The row order in the names files is identical 
to the column order in the kraken2 table. 

Ancient dental calculus: `anc_calc.kraken2.gtdb.report_mpa.tsv.gz`  
Baboon calculus: `baboon_calc.kraken2.gtdb.report_mpa.tsv`; `baboon_names.tsv`  
Chimpanzee calculus: `chimp_calc.kraken2.gtdb.report_mpa.tsv`  
Gorilla calculus: `gorilla_calc.kraken2.gtdb.report_mpa.tsv`; `gorilla_names.tsv`  
Howler monkey calculus: `howler_calc.kraken2.gtdb.report_mpa.tsv`; `howler_names.tsv`  
Industrial buccal mucosa: `ind_buccal_mucosa.kraken2.gtdb.report_mpa.tsv`; `indbm_names.tsv`  
Industrial dental plaque: `ind_plaque.kraken2.gtdb.report_mpa.tsv`; `indpq_names.tsv`  
Industrial saliva: `ind_saliva.kraken2.gtdb.report_mpa.tsv`; `indsal_names.tsv`  
Modern human dental calculus: `mod_calc.kraken2.gtdb.report_mpa.tsv`; `mod_calc_names.tsv`  
Non-human primate oral swabs: `nhp_os.kraken2.gtdb.report_mpa.tsv`; `nhp_os_names.tsv`  
Non-industrial buccal mucosa: `nonind_buccal_mucosa.kraken2.gtdb.report_mpa.tsv`; `nonindbm_names.tsv`  
Non-industrial dental plaque: `nonind_plaque.kraken2.gtdb.report_mpa.tsv`; `nonindpq_names.tsv`  
Non-industrial saliva: `nonind_saliva.kraken2.gtdb.report_mpa.tsv`; `nonindsal_names.tsv`  
Non-oral samples for preservation assessment: `sourceTracker_sources.kraken2.gtdb.report_mpa.tsv`; `source_names.tsv`  

Tables of genes enriched in different sample types are in `05-results`.  
Human non-shedding surface vs. human shedding surface: `pq_calc_vs_bms_sig_genes.tsv`  
Human shedding surface vs. non-human primate shedding surface: `bms_nhpos_sig_genes.tsv`  
Human non-shedding surface vs. non-human primate non-shedding surface: `calc_h_vs_nhp_sig_genes.tsv`  
Human ancient non-shedding vs. human modern non-shedding: `calc_anc_vs_mod_sig_genes.tsv.gz`  

