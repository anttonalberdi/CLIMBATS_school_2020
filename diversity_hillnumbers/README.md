# Diversity analyses using Hill numbers
Hill numbers provide a general statistical framework that is sufficiently robust and flexible to address a wide range of scientific questions that molecular ecologists regularly try to answer through measurement, estimation, partitioning, and comparison of diversities (Jost 2006; Tuomisto 2010a; Chao et al. 2014a). This lesson aims at introducing the basics of Hill numbers to perform different types of diversity analyses from OTU or ASV tables generated from molecularly characterised biological systems: e.g. diet, microbiomes.

This lessons is largely based on an article published in Molecular Ecology Resources, in which we provided  guidelines for implementing Hill numbers in diversity analyses of metabarcoding data. 

Alberdi A, Gilbert MTP. (2019). A guide to the application of Hill numbers to DNA‐based diversity analyses. *Molecular Ecology Resources*. 19(4): 804-817. [https://doi.org/10.1111/1755-0998.13014](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13014)

## Table of contents

1. [Diversity metrics](1-Diversity_metrics.md)
2. [Introduction to Hill numbers](2-Introduction_to_Hill_numbers.md)
3. [Understanding Hill numbers](3-Understanding_Hill_numbers.md)
4. [Diversity partitioning using Hill numbers](4-Diversity_partitioning_and_dissimilarity_measurement.md)
5. [Integral diversity analysis using a real dataset](5-Integral_diversity_analysis_using_a_real_dataset.md)

## The package hilldiv

This lesson relies largely on the R package [hilldiv](https://github.com/anttonalberdi/hilldiv). hilldiv is an R package that provides a set of functions to assist analysis of diversity for diet reconstruction, microbial community profiling or more general ecosystem characterisation analyses based on Hill numbers, using OTU/ASV tables and associated phylogenetic trees as inputs. The package includes functions for (phylo)diversity measurement, (phylo)diversity profile plotting, (phylo)diversity comparison between samples and groups, (phylo)diversity partitioning and (dis)similarity measurement. All of these grounded in abundance-based and incidence-based Hill numbers.

### Installation

To install **hilldiv** in your R environment you can rely on the built-in install.packages() function:
````R
install.packages("hilldiv")
library(hilldiv)
library(geiger)
library(ape)
````
If you want to use the latest development version available at Github, you need to 1) install devtools, 2) load devtools library, 3) install **hilldiv** using devtools and 4) finally load **hilldiv** library to your environment.

````R
install.packages("devtools")
library(devtools)
remove.packages("hilldiv") #if already an older version installed
install_github("anttonalberdi/hilldiv")
library(hilldiv,quietly=TRUE)
````
If not installed, it will automatically install the following dependencies: ggplot2, ggpubr, RColorBrewer, data.table, ape, vegan, geiger, qgraph and FSA.

If the console returns the following error:
````R
"tar: Failed to set default locale"
````
Type the following in the console and restart R.
````R
system('defaults write org.R-project.R force.LANG en_US.UTF-8')
````
