# treeClust-benchmarking
R scripts and data for treeClust benchmarking manuscript

This is a brief explanation of the R scripts used for the benchmarking of the treeClust algorithm. This work is documented in more detail in our manuscript (Kustatscher, Grabowski and Rappsilber, https://www.biorxiv.org/content/10.1101/578971v1). The treeClust algorithm was developed by Buttrey & Whitaker ("treeClust: an R package for tree-based clustering dissimilarities", The R Journal, 2015). All scripts were written in Ubuntu 18.04 using R version 3.5.1 ("Feather Spray").

#### treeClust_bench_relationship_types.R
This script compares treeClust and correlation measures for their ability to detect non-linear relationships. It produces figure 1 of the manuscript.

#### treeClust_benchmarking.R
This script creates a series of synthetic datasets with a range of different properties and tests how treeClust and correlation measures respond to that. It produces figure 2 of the manuscript. Note that execution of this script will take a few hours on a standard PC.

#### Benchmark_analysis_ProHD.R
This script tests if our findings on synthetic data hold true for ProteomeHD, a real proteomics dataset that we published previously {UPDATE REF}. In addition, this script uses the file Reactome_TP_FP.csv, which contains a set of true and false positive protein pairs which we designed previously based on the Reactome database {UPDATE REF}. Essentially, true positive interactions are proteins that function in the same detailed pathway (excluding pathways with more than 200 members), while false positives are proteins that were annotated by Reactome but not to the same pathway. Both files (ProteomeHD_v1_1.csv and Reactome_TP_FP.csv) are in this repository in compressed format. To run this script, extract the two input files and make sure R's working directory is set to the folder containing the files. Running the script will take several hours on a standard PC. It produces figures 3, 4 and S1.   
