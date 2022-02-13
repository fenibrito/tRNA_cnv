# tRNA_cnv
 
README

This repository contains all the data and code used in the project

The `tRNA_cnv` directory holds all the intermediate files (`data`) and the R scripts (`scripts`) needed to generate all the Figures and Tables in the paper:


`/home/feni/repository/trna_project/scripts/download_filter_gtrnadb_data.R` contains the code to download the GtRNAdb data and to make all the modifications needed for further analysis
`/home/feni/repository/trna_project/scripts/taxID.R` contains the code and instructions to recover the taxonomic lineage for each specie on the GtRNAdb

**running the first two scripts is necessary for the proper execution of the following scripts**


`/home/feni/repository/trna_project/scripts/trna_count_per_domain.R` contains the code to reproduce Figures 1A and 1B
`/home/feni/repository/trna_project/scripts/anticodon_per_genome_domain.R` contains the code to reproduce Figures 1C
`/home/feni/repository/trna_project/scripts/anticodon_per_genome_taxa.R` contains the code to reproduce Figures 1D
`/home/feni/repository/trna_project/scripts/mean_trna_isotype_domain.R` contains the code to reproduce Figures 1E
`/home/feni/repository/trna_project/scripts/mean_trna_anticodon_domain.R` contains the code to reproduce Figures 1F
`/home/feni/repository/trna_project/scripts/codon_usage.R` contains the code to reproduce Figures 2A and 2B
`/home/feni/repository/trna_project/scripts/anticodon_obs_freq.R` contains the code to reproduce Figures 2C, 2D and 2E
`/home/feni/repository/trna_project/scripts/genome_size_vs_gene_count.R` contains the code to reproduce Figures 2F, 2G, 2H, 3A, 3B, 3C, 3D and 3E
`/home/feni/repository/trna_project/scripts/n_anticodon_vs_genomic_features.R` contains the code to reproduce Figure 3F, 3G and H
`/home/feni/repository/trna_project/scripts/compare_genome_information.R` contains the code to reproduce Figure 3I
`/home/feni/repository/trna_project/scripts/trna_gene_length.R` contains the code to reproduce Figures 3J
`/home/feni/repository/trna_project/scripts/intron_count.R` contains the code to reproduce Figures 3J
`/home/feni/repository/trna_project/scripts/cnv_interspecie.R` contains the code to reproduce Figure 4A and 4B
`/home/feni/repository/trna_project/scripts/trfs_per_sp.R` contains the code to reproduce Figure 5A,5B,5C,5D,5F,5G,5H and 5I
`/home/feni/repository/trna_project/scripts/scatterplot3d_isoaccep_anticodon_mean_dom.R` contains the code to reproduce Supplementary Figure 1