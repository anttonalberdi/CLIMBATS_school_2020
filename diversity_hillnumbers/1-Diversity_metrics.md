# Diversity metrics

## Defining types for diversity quantification
In community ecology, individuals (i.e., recorded entities) have been traditionally classified into taxonomic species (i.e., types). Therefore, diversity measurements have commonly been carried out at species level (e.g., species richness and species diversity), principally as determined based on morphological features (MacArthur, 1965; Pielou, 1966). The implementation of DNA‐based molecular approaches now enables (in principle) diversity to be measured at a much finer scale—that of DNA sequence variation. Molecularly defined types, broadly known as **OTUs or MOTUs** (molecular operational taxonomic units, Blaxter et al., 2005), are becoming the preferred types with which to quantify diversity in many fields of the biological sciences. When using molecular approaches, the **recorded entities are no longer individuals, but DNA sequences**, and the classification into types is not any longer based on morphological features, but the level of dissimilarity between DNA sequences.

![Traditional and molecular definition of types](https://github.com/anttonalberdi/CLIMBATS_school_2020/blob/master/diversity_hillnumbers/images/types.png)

## Abundance-based vs. incidence-based quantification of diversity
Diversity measurements require assignment of an importance value to each of the detected types. In traditional community ecology, this has been done using metrics such as individual counts, biomass or spatial units, depending on the type of system, research question and fieldwork strategy. Molecular analyses provide a different type of data that could provide such information, namely the amount of DNA sequences assigned to each OTU (Deagle et al., 2019). 

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

## Diversity components and traditional diversity metrics 
Biological diversity is a complex feature that can be decomposed into richness, evenness and regularity components (Alberdi et al. 2019; Jost 2010). Each of these components measure different properties of the diversity, and can be shaped by different ecological forces (Wilsey and Stirling 2007).

### Richness
The simplest measure of diversity is OTU richness (McIntosh, 1967). As this only considers whether an OTU is present or absent in the system, abundant and rare OTUs are given the same weight. 

````R
library(hilldiv)

# Define an even system
evensystem <- c(5,5)
names(evensystem) <- c("OTU1","OTU2")
evensystem
# OTU1 OTU2 
#    5    5 

#Compute richness
index_div(evensystem,index="richness")
# 2

# Define an uneven system
unevensystem <- c(999,1)
names(unevensystem) <- c("OTU1","OTU2")
unevensystem
# OTU1 OTU2 
#  999    1 

#Compute richness
index_div(example1,index="richness")
# 2
````
However, the multiple OTUs present in a system are seldom distributed evenly; thus, richness is rarely the best approach with which to reflect the diversity of a system. Consider for instance, a simple system characterized with 1,000 sequence reads, in which 990 belong to OTU1 and 10 to OTU2. This would yield a richness value of 2, even though the system is overwhelmingly dominated by OTU1. 

### Evenness: Shannon and Simpson indices

However, in real biological systems OTUs are rarely evenly represented, thus including the evenness component, which measures the balance of the relative representation of each OTU, yields more meaningful results. Two popular metrics that that account for both OTU richness and evenness are the Shannon and Simpson indices. These two metrics differ in the way in which the importance of abundant and rare OTUs are weighed. While the former weights each OTU for its quantitative representation —either number DNA sequences or samples, the latter overweights the importance of abundant (OTUs (Alberdi and Gilbert 2019).  It is critical to note, however, that unlike richness, neither the Shannon index nor the Simpson index are actual measures of diversity. The former measures entropy thus yields the uncertainty in the OTU identity of a randomly chosen sequence in the system. The latter provides the probability that two randomly chosen DNA sequences actually belong to different OTUs (Chao, Chiu, et al., 2014a). Consequently, the values that Shannon and Simpson indices yield are difficult to interpret—the values in the previous example are 0.079 and 0.020, respectively, and do not exhibit the intuitive properties ecologists expect from a diversity measurement.

````R
library(hilldiv)

unevensystem <- c(999,1)
names(unevensystem) <- c("OTU1","OTU2")

index_div(unevensystem,index="shannon")
# 0.007907255
index_div(unevensystem,index="simpson")
# 0.001998
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

### Regularity: Faith’s PD,  Allen’s H and Rao's Q

A third component with key biological relevance, yet often overlooked in metabarcoding, is regularity, i.e. the degree of similarity across the OTUs. Richness, Shannon index and Simpson index treat OTUs as independent elements, thus overlooking that in most biological systems OTUs tend to be functionally, phylogenetically or ecologically correlated. Faith’s PD,  Allen’s H and Rao's Q are three metrics that incorporate the regularity component, the first considering only richness and regularity, while the last two accounting for  richness, evenness and regularity. 

The relationship between traditional diversity indices and phylogenetic indices can be observed when comparing a star-like even tree in which the distance between all OTUs is identical with an uneven tree in which the distances between OTUs differ.  

````R
#Create an even system
evensystem <- c(5,5,5)
names(evensystem) <- c("OTU1","OTU2","OTU3")

#Create trees
library(phytools)
eventree <- starTree(c("OTU1","OTU2","OTU3"), branch.lengths=c(1,1,1))
uneventree <- read.tree(text="((OTU1:0.5,OTU2:0.5):0.5,OTU3:1);")

#Plot trees
plot(eventree)
plot(uneventree)
````

Faith's PD is identical to richness when the tree that relates all three OTUs is star-like. However the phylogenetic richness value decreases when using the unven tree.

````R
#Compute Faith's PD
index_div(evensystem, index="richness")
# [1] 3
index_div(evensystem, tree=eventree, index="faith")
# [1] 3
index_div(evensystem, tree=uneventree, index="faith")
# [1] 2.5
````

Allen's H is identical to Shannon index when the tree that relates all three OTUs is star-like. However the phylogenetic diversity value decreases when using the unven tree.

````R
#Compute Allen's H
index_div(evensystem, index="shannon")
# [1] 1.098612
index_div(evensystem, tree=eventree, index="allen")
# [1] 1.098612
index_div(evensystem, tree=uneventree, index="allen")
# [1] 0.8675632
````

Rao's Q is identical to Simpson index when the tree that relates all three OTUs is star-like. However the phylogenetic diversity value decreases when using the unven tree.

````R
#Compute Rao's Q
index_div(evensystem, index="simpson")
# [1] 0.6666667
index_div(evensystem, tree=eventree, index="rao")
# [1] 0.6666667
index_div(evensystem, tree=uneventree, index="rao")
# [1] 0.5555556
````
