# Diversity analyses using Hill numbers
Hill numbers provide a general statistical framework that is sufficiently robust and flexible to address a wide range of scientific questions that molecular ecologists regularly try to answer through measurement, estimation, partitioning, and comparison of diversities (Jost 2006; Tuomisto 2010a; Chao et al. 2014a). In an article published in Molecular Ecology Resources, we provided some guidelines for implementing Hill numbers in diversity analyses of metabarcoding data. The aim of this wiki is to serve as a continuation of that article with a more practical focus.

Alberdi A, Gilbert MTP. (2019). A guide to the application of Hill numbers to DNA‐based diversity analyses. *Molecular Ecology Resources*. 19(4): 804-817. [https://doi.org/10.1111/1755-0998.13014](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13014)

All the examples shown in this tutorial rely on the R package hilldiv, which can be downloaded from CRAN.

````R
install.packages("hilldiv")
library(hilldiv)
````

## Table of contents

1. Introduction
2. [Diversity metrics](2-Diversity_metrics.md)
3. Hill numbers

