# Diversity analyses using Hill numbers
Hill numbers provide a general statistical framework that is sufficiently robust and flexible to address a wide range of scientific questions that molecular ecologists regularly try to answer through measurement, estimation, partitioning, and comparison of diversities (Jost 2006; Tuomisto 2010a; Chao et al. 2014a). In an article published in Molecular Ecology Resources, we provided some guidelines for implementing Hill numbers in diversity analyses of metabarcoding data. The aim of this wiki is to serve as a continuation of that article with a more practical focus.

Alberdi A, Gilbert MTP. (2019). A guide to the application of Hill numbers to DNA‐based diversity analyses. *Molecular Ecology Resources*. 19(4): 804-817. [https://doi.org/10.1111/1755-0998.13014](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.13014)

All the examples shown in this tutorial rely on the R package hilldiv, which can be downloaded from CRAN.

````R
install.packages("hilldiv")
library(hilldiv)
````

## Defining types for diversity quantification
In community ecology, individuals (i.e., recorded entities) have been traditionally classified into taxonomic species (i.e., types). Therefore, diversity measurements have commonly been carried out at spe‐ cies level (e.g., species richness and species diversity), principally as determined based on morphological features (MacArthur, 1965; Pielou, 1966). The implementation of DNA‐based molecular approaches now enables (in principle) diversity to be measured at a much finer scale—that of DNA sequence variation. Molecularly defined types, broadly known as OTUs or MOTUs (molecular operational taxonomic units, Blaxter et al., 2005), are becoming the preferred types with which to quantify diversity in many fields of the biological sciences. When using molecular approaches, the recorded entities are no longer individuals, but DNA sequences, and the classification into types is not any longer based on morphological features, but the level of dissimilarity between DNA sequences.

Diversity measurements require assignment of an importance value to each of the detected types. In traditional community ecology, this has been done using metrics such as individual counts, biomass or spatial units, depending on the type of system, research question and fieldwork strategy. Molecular analyses provide a different type of data that could provide such information, namely the amount of DNA sequences assigned to each OTU (Deagle et al., 2019). 

## Abundance-based vs. incidence-based quantification of diversity
The diversity of molecularly characterised biological systems can be measured using two main approaches: abundance-based and incidence-based. Although incidence data are less informative than abundance data, it is both easier to collect, more comparable, and has been extensively used under the niche theory framework. When dealing with DNA‐derived data, incidence‐based approaches have particular relevance, given the limited quantitative relationship that exists between the biomass in the actual system and the DNA sequences produced (Lamb et al., 2019), which might challenge the representativeness of abundance data. However, consensus has not been reached within the molecular ecology research community about which approach is the most appropriate, as simulations have shown that analyses based on incidence data often overestimate the importance of rare OTUs, and abundance data might provide a more accurate view of the diversity even with moderate recovery biases (Deagle et al., 2019). In fact, we recently showed that, when sample sizes and the effect sizes are large enough, both approaches yield biologically meaningful results (Alberdi et al. 2020).

While either approach might be valid depending on the research question and the study design, and traditional diversity indices as well as Hill numbers can be computed on both types of data, it is important to acknowledge the basic differences between abundance‐based and incidence‐based diversity metrics. 

In the **abundance‐based** approach, the unit used to compute diversity is the **count of DNA sequences** assigned to each OTU. 

In the **incidence‐based** approach, the unit used to compute diversity is the **count of samples** in which an OTU is present. 

### Transforming abundance into incidence data
Abundance data (e.g. OTU table) can be easily transformed into incidence data using the hilldiv function to.incidence().

````R
library(hilldiv)

# Create an abundance table, with columns defining samples, and rows OTUs. 
abundance.table <- cbind(Sample1=c(8,2,20),Sample2=c(4,0,10),Sample3=c(0,0,12),Sample4=c(2,6,0))
rownames(abundance.table) <- c("OTU1","OTU2","OTU3")
abundance.table 
#      Sample1 Sample2 Sample3 Sample4
# OTU1       8       4       0       2
# OTU2       2       0       0       6
# OTU3      20      10      12       0

# Transform the abundance table into incidence
incidence.table <- to.incidence(abundance.table)
incidence.table
# OTU1 OTU2 OTU3 
#    3    2    3 
# The bi-dimensional abundance data composed of 4 samples has been collapsed into a 
# uni-dimensional indicence vector, showing that OTU1 is present in 3 samples, OTU2 
# in 2 samples and OTU3 in 3 samples.
````

The transformation from abundance to incidence data is usually performed by groups. For example, if the aim is to compare the diet of species, abundance data of individual samples can be collapsed into incindence data per species. This can be achieved by specifying a so-called hierarchy table that shows the relationship between samples and groups. As many groups as wished can be defined, but each sample can only be related to one of the groups.

````R
library(hilldiv)
#Define a hierarchy table grouping samples into two species
hierarchy.table <- cbind(Samples=c("Sample1","Sample2","Sample3","Sample4"),Groups=c("Species1","Species1","Species2","Species2"))
hierarchy.table
#      Samples   Groups    
# [1,] "Sample1" "Species1"
# [2,] "Sample2" "Species1"
# [3,] "Sample3" "Species2"
# [4,] "Sample4" "Species2"

#Transform the abundance table into an incidence table considering the grouping information
incidence.table <- to.incidence(abundance.table, hierarchy=hierarchy.table)
incidence.table
#     Species1 Species2
# OTU1        2        1
# OTU2        1        1
# OTU3        2        1
````

## Diversity components
Biological diversity is a complex feature that can be decomposed into richness, evenness and regularity components (Alberdi et al. 2019; Jost 2010). Each of these components measure different properties of the diversity, and can be shaped by different ecological forces (Wilsey and Stirling 2007).

### Richness
The simplest measure of diversity is OTU richness (McIntosh, 1967). As this only considers whether an OTU is present or absent in the system, abundant and rare OTUs are given the same weight. However, the multiple OTUs present in a system are seldom distributed evenly; thus, richness is rarely the best approach with which to reflect the diversity of a system. Consider for instance, a simple system characterized with 1,000 sequence reads, in which 990 belong to OTU1 and 10 to OTU2. This would yield a richness value of 2, even though the system is overwhelmingly dominated by OTU1. 

````R
library(hilldiv)

example1 <- c(999,1)
index_div(example1,index="richness")
2
````
### Traditional diversity indices and their limitations
Thus, metrics such as the Shannon or the Simpson indices, which also account for the evenness of the system, are con‐ sidered more representative of the diversity of a system. It is critical to note, however, that unlike richness, neither the Shannon index nor the Simpson index are actual measures of diversity. The former measures entropy thus yields the uncertainty in the OTU identity of a randomly chosen sequence in the system. The latter provides the probability that two randomly chosen DNA sequences actually belong to different OTUs (Chao, Chiu, et al., 2014a). Consequently, the values that Shannon and Simpson indices yield are difficult to interpret—the values in the previous example are 0.079 and 0.020, respectively, and do not exhibit the intuitive properties ecologists expect from a diversity measurement.


````R
library(hilldiv)

example1 <- c(999,1)
index_div(example1,index="shannon")
0.007907255
index_div(example1,index="simpson")
0.001998
````

Specifically, our intuitive notion of diversity would expect that when doubling the number of OTUs in a system, then the diversity measured should also double. This is known as the “doubling property” or “replication principle” (Chao, Chiu, & Jost, 2010; Hill, 1973; Jost, 2006). For example, if the diet of one bat species is comprised of 15 moth species with even abundances, and the diet of another species encompasses 30 moths also with even abundances, intuitively we would conclude that the second bat's diet is twice as diverse (100% more diverse) as the first one. 

While richness owns that property, most diversity indices do not. The Shannon entropy only increases from 2.7 (15 species) to 3.4 (30 species), which might suggest a diversity gain of 26%, and the Simpson index only increases from 0.93 to 0.96, which might suggest a gain of barely 3%. Hence, treating diversity indices as diversity values has noticeable practical consequences, as they all vary in range and behaviour (Jost, 2006).

````R
library(hilldiv)

bat1 <- c(rep(1,15))
bat1
 [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
bat2 <- c(rep(1,30))
bat2
 [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 
#Richness
index_div(bat1,index="richness")
15
index_div(bat2,index="richness")
30

#Shannon
index_div(bat1,index="shannon")
2.70805
index_div(bat2,index="shannon")
3.401197

#Simpson
index_div(bat1,index="simpson")
0.9333333
index_div(bat2,index="simpson")
0.9666667
````

### Hill numbers
Richness, Shannon index and Simpson index belong to a single statistical framework, as they all are monotonic functions of the basic sum q𝜆=ΣSi=1pqi, that is, the sum of the relative abundances of the types (pi) elevated to the q value (Jost, 2006; Keylock, 2005). 

This implies that Hill numbers (qD), or actual diversities, rather than entropies (e.g. Shannon index) or probabilities (e.g. Simpson index), can be formulated in terms of the basic sum (qλ) and the parameter q.

**Hill number expression**

This expression was first discovered by Hill (1973), hence the use of the name “Hill numbers” to refer to the output of this formula. Hill numbers have two major advantages over diversity indi‐ ces: (a) the interpretation of the measure and its measurement unit is always the same (Chao, Chiu, et al., 2014a; Tuomisto, 2010a), and ii) the sensitivity towards abundant and rare OTUs can be modulated with the parameter q. 

#### Effective number of OTUs
The expression yields a diver‐ sity measure in “effective number of OTUs”, that is, the number of equally abundant OTUs that would be needed to give the same value of diversity (Hill, 1973; Jost, 2006). When all OTUs in a system have the same relative abundances, as in the moth example given above, the effective number of OTUs for all q values equals the actual number of OTUs, namely richness.

````R
library(hilldiv)

# Create an even system composed of 30 OTUs each represented by 1 sequence
evensystem <- c(rep(1,30))
evensystem
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

# When the types are evenly distributed in the system, Hill numbers 
# of q-value 0, 1 and 2 (or any other q value) equal richness

index_div(evensystem,index="richness")
# [1] 30
hill_div(evensystem,qvalue=0)
# [1] 30
hill_div(evensystem,qvalue=1)
# [1] 30
hill_div(evensystem,qvalue=2)
# [1] 30
````

When the relative abundances of the types vary however, then the effective number of OTUs for q > 0 values decreases progresivelly. The higher the heterogeneity between types, the sharper will be decrease of the effective number of OTUs.

````R
library(hilldiv)

# Create an uneven system composed of 15 OTUs represented 
# by 1 sequence and 15 OTUs represented by 5 sequences
unevensystem <- c(rep(1,15),rep(5,15))
unevensystem
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5

# When types (OTUs) are not evenly distributed in the system, Hill numbers 
# of q-value > 0 show a decreasing effective number of types (OTUs)

index_div(unevensystem,index="richness")
# [1] 30
hill_div(unevensystem,qvalue=0)
# [1] 30
hill_div(unevensystem,qvalue=1)
# [1] 23.53793
hill_div(unevensystem,qvalue=2)
# [1] 20.76923
hill_div(unevensystem,qvalue=5)
# [1] 18.83793
````
In extreme cases in which the system is dominated by a few equally abundant OTUs, the effective number of OTUs will approach the number of those abundant OTUs.

````R
library(hilldiv)

# Create a super uneven system composed of 29 OTUs represented 
# by 1 sequence and 1 OTU represented by 971 sequences
superunevensystem <- c(rep(1,29),971)
superunevensystem
# [1]   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
# [25]   1   1   1   1   1 971

# When types (OTUs) are very unevenly distributed in the system, Hill numbers 
# of q-value > 0 show a sharp drop of effective number of types (OTUs),
# approaching the number of dominant types (OTUs)

index_div(superunevensystem,index="richness")
# [1] 30
hill_div(superunevensystem,qvalue=0)
# [1] 30
hill_div(superunevensystem,qvalue=1)
# [1] 1.257225
hill_div(superunevensystem,qvalue=2)
# [1] 1.060592
hill_div(superunevensystem,qvalue=5)
# [1] 1.037471
````
#### Plotting diversity profiles
Hill numbers also enable diversity profiles of systems and subsystems to be plotted as continuous functions of the parameter q. This is useful to characterize the OTU abundance distri‐ bution of a system, as different compositions and abundance distributions can yield the same value for a particular order of diversity (e.g., q=1),but not for many of them(e.g. q=0, q=0.5 and q=1). Hill numbers convey all information contained in a species abundance distribution at a glance (Chao, Chiu, et al., 2014a; Leinster & Cobbold, 2012).

````R
library(hilldiv)

#Define the three model systems 
evensystem <- c(rep(1,30))
unevensystem <- c(rep(1,15),rep(5,15))
superunevensystem <- c(rep(1,29),971)

#Merge them in a single OTU table

mergedsystem <- cbind(evensystem,unevensystem,superunevensystem)
mergedsystem_profile <- div_profile(mergedsystem)
div_profile_plot(mergedsystem_profile)
````

#### Modulating the q value to adjust to study systems' properties
As shown above, the sensitivity towards abundant and rare OTUs can be modulated using the scaling parameter q, known as the “order” of diversity (Jost, 2006). **The larger the q value, the higher the importance attributed to abundant OTUs.** The ability to modulate the sensitivity towards abundant and rare OTUs by modifying a single parameter (q) is a useful means with which to adjust diversity measurements to the type of data and re‐ search question. For example, when rare types are considered to be of low importance (e.g., when attempting to define a core diet or mi‐ crobiome), or when rare types are considered untrustworthy due to technical issues (e.g., PCR or sequencing errors), researchers might opt for using a high q value, for example, q = 2, which overweighs abundant OTUs. The result can be interpreted as the effective num‐ ber of dominant OTUs in the system (Chao, Chiu, et al., 2014a). In contrast, if rare types are considered essential for the system, or re‐ searchers do not trust the relative abundance data due to potential technical biases, researchers might opt for using a q value of 0 that simply counts the number of types.

#### Relationship between Hill numbers and traditional popular indices
Three q values are particularly relevant, both for their significance, and their close relationship to popular diversity indices: q = 0, q = 1 and q = 2. Indeed, common diversity indices can be transformed to Hill numbers by applying simple mathematical transformations, as shown below.

As we have observed in the previous example, when a diversity of order zero (q = 0) is applied to the Hill numbers expression, it becomes insensitive to OTU frequencies, thus yielding a richness value. As the relative abundances of OTUs are overlooked, rare OTUs are overweighed.

````R
index_div(superunevensystem,index="richness")
# [1] 30
hill_div(superunevensystem,qvalue=0)
# [1] 30
````

A q value of 1 (in practical terms its limit, as the Hill number is undefined for q = 1) is the value that weighs OTUs by their frequency, without disproportionately favouring either rare or abundant ones (Jost, 2006). The value it yields is exactly the exponential of the Shannon index. In fact, q values under unity favour rare OTUs, while values above one favour abundant OTUs (Keylock, 2005). 

````R
hill_div(superunevensystem,qvalue=1)
# [1] 1.257225
index_div(superunevensystem,index="shannon")
# [1] 0.2289003
exp(index_div(superunevensystem,index="shannon"))
# [1] 1.257217
````

When a q value of 2 is applied, abundant OTUs are overweighed, and the formula yields the multiplicative inverse of the Simpson index.

````R
hill_div(superunevensystem,qvalue=2)
# [1] 1.060592
index_div(superunevensystem,index="simpson")
# [1] 0.05713
1/(1-index_div(superunevensystem,index="simpson"))
# [1] 1.060592
````

