# HOCI_SRT_R_code
Code to run calculations for the SARS-CoV-2 sequence reporting tool as described in https://www.medrxiv.org/content/10.1101/2020.11.12.20230326v1

This repository includes:
- Dummy_data_SRT_R_script.R: R script to read in data and run the analysis
- geo_cluster_fun_SNP_para.R: R function for geographic clustering model (can be run with multiple cores if needed, as this can be slow for large datasets)
- gen_posterior_fun_SNP_M2_locmatch_2.R: R function to carry out SRT calculations
- hociSiteMetadata_dummy.csv: Artificially constructed and named meta-data, with variables named and formatted to match the COV-GLUE HOCI tool
- hociSiteSequences_MSA_dummy: Sequences to match the dummy meta-data (these need to be aligned, unlike the COV-GLUE HOCI tool)
- dataset_focus_results.csv: Output of the R scripts, with SRT results reported for focus sequences
