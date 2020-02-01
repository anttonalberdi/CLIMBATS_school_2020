# Diversity analyses using Hill numbers
Texts have been modified from:

Alberdi A, Gilbert MTP. 2019. A guide to the application of Hill numbers to DNA based diversity analyses. Molecular Ecology Resources 19: 804-817.

Researchers often need to quantify how diverse different systems are, for example, to assess ecosystem functioning (Cardinale, Palmer, & Collins, 2002) or to measure any species' niche breadth (Forister et al., 2015). It is also common to compare the composition of different systems, under experimental set‐ups to measure the dif‐ ferences yielded by different treatments (e.g., Gevers et al., 2014), or in observational designs to assess whether and how much dietary niches differ (e.g., Kartzinel et al., 2015). A myriad of approaches and tools has been developed over the last century to perform such operations, each embedded within a statistical background, with popular examples including Shannon index quadratic entropy (Rao, 1982), Pianka's niche overlap (Pianka, 1974) and Unifrac distances (Lozupone & Knight, 2005).

Regardless of the specific approach implemented, researchers need to make at least four essential choices when analysing the diversity of a biological system. 

* First, they must define the unit that encompasses biologically alike entities to be used to measure diversity, hereafter referred to as the “type”. Although community ecologists have traditionally measured diversity using the taxonomic species as the type (Pielou, 1966), with the implementation of molecular approaches, this is no longer a general rule (Blaxter et al., 2005). 

* Second, it is necessary to consider how detections of these types are treated, either as presence/absence (=incidence) or quantitatively (=abundance), and if the latter, how abundant and rare detections are weighed (Jost, 2006). 

* Third, researchers need to con‐ sider whether, and how, the phylogenetic—or ecological—relations between detected organisms will be accounted for when measur‐ ing diversity (Chao, Chiu, & Jost, 2014a). 

* Finally, researchers should ideally assess whether the data are representative of the biological system studied, and if needed, take the necessary measures to cor‐ rect the biases.

## Defining types for diversity quantification

In community ecology, individuals (i.e., recorded entities) have been traditionally classified into taxonomic species (i.e., types). Therefore, diversity measurements have commonly been carried out at spe‐ cies level (e.g., species richness and species diversity), principally as determined based on morphological features (MacArthur, 1965; Pielou, 1966). The implementation of DNA‐based molecular approaches now enables (in principle) diversity to be measured at a much finer scale—that of DNA sequence variation. Molecularly defined types, broadly known as OTUs or MOTUs (molecular operational taxonomic units, Blaxter et al., 2005), are becoming the preferred types with which to quantify diversity in many fields of the biological sciences. When using molecular approaches, the recorded entities are no longer individuals, but DNA sequences, and the classification into types is not any longer based on morphological features, but the level of dissimilarity between DNA sequences.

## Weighing the importance of types
Diversity measurements require assignment of an importance value to each of the detected types. In traditional community ecology, this has been done using metrics such as individual counts, biomass or spatial units, depending on the type of system, research question and fieldwork strategy. Molecular analyses provide a different type of data that could provide such information, namely the amount of DNA sequences assigned to each OTU (Deagle et al., 2019). There are multiple approaches that enable differential weighing of abundant and rare OTUs.

### Richness
The simplest measure of diversity is OTU richness (McIntosh, 1967). As this only considers whether an OTU is present or absent in the system, abundant and rare OTUs are given the same weight. However, the multiple OTUs present in a system are seldom distributed evenly; thus, richness is rarely the best approach with which to reflect the diversity of a system. Consider for instance, a simple system characterized with 1,000 sequence reads, in which 990 belong to OTU1 and 10 to OTU2. This would yield a richness value of 2, even though the system is overwhelmingly dominated by OTU1. 

````R
library(hilldiv)

example1 <- c(999,1)
index_div(example1,index="richness")
2
````
### Eveness

Thus, metrics such as the Shannon or the Simpson indices, which also account for the evenness of the system, are con‐ sidered more representative of the diversity of a system. It is critical to note, however, that unlike richness, neither the Shannon index nor the Simpson index are actual measures of diversity. The former measures entropy thus yields the uncertainty in the OTU identity of a randomly chosen sequence in the system. The latter provides the probability that two randomly chosen DNA sequences actually belong to different OTUs (Chao, Chiu, et al., 2014a). Consequently, the values that Shannon and Simpson indices yield are difficult to interpret—the values in the previous example are 0.079 and 0.020, respectively, and do not exhibit the intuitive properties ecologists expect from a diversity measurement.


````R
library(hilldiv)

example1 <- c(999,1)
index_div(example1,index="shannon")
0.007907255
index_div(example1,index="simpson")
0.001998
````

