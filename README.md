# treeClust-benchmarking
R scripts for treeClust benchmarking manuscript

This is a brief explanation of the R scripts used for the benchmarking of the treeClust algorithm. This work is documented in more detail in our manuscript (ADD bioRxiv LINK HERE). The treeClust algorithm was developed by Buttrey & Whitaker ("treeClust: an R package for tree-based clustering dissimilarities", The R Journal, 2015).


#### treeClust_bench_relationship_types.R
This script compares treeClust and correlation measures for their ability to detect non-linear relationships. It produces figure 1 of the manuscript.

#### treeClust_benchmarking.R
This script creates a series of synthetic datasets with a range of different properties and tests how treeClust and correlation measures respond to that. It produces figure 2 of the manuscript. Note that execution of this script will take a few hours on a standard PC.

