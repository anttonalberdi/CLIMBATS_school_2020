# Integral diversity analysis using a real dataset
This lesson aims at showcasing the use of Hill numbers to analyse multiple aspects of the diversity of a real dataset. The data is a subset of the bat dietary data generated under the following publication:

Alberdi A, Razgour O, Aizpurua O, Novella-Fernandez R, Aihartza J, Budinski I, Garin I, Ibáñez C, Izagirre E, Rebelo H, Russo D, Vlaschenko A, Zhelyazkova V, Zrncic V, Gilbert MTP. 2020. **DNA metabarcoding and spatial modelling link diet diversification with distribution homogeneity in European bats.** Nature Communications. Accepted.

The dataset has been downscaled and slightly modified for educational purposes. The raw material used to generate the data can be retrieved from the Zenodo repository: [https://zenodo.org/record/3610756#.XjhGGxNKiCM](https://zenodo.org/record/3610756#.XjhGGxNKiCM)

## Setting the working environment
First of all download the files stored in the [data directory](https://github.com/anttonalberdi/CLIMBATS_school_2020/tree/master/diversity_hillnumbers/data) of this lesson, and save them in the local directory you will use for the analyses.

````R
#The working directory in this example is my Desktop
setwd("/Users/anttonalberdi/Desktop")
````

Load all the libraries we need for performing the analyses.
````R
library(hilldiv)
library(ape)
````

Load the data files that will be used in the lesson
````R
#Low the OTU table
otutable <- read.csv("EuropeBatDiet.csv",row.names=1)
#Visualise the first 6 rows
head(otutable)

#Low the OTU table
tree <- read.tree("EuropeBatDiet.tree")
#Visualise it
plot(tree)

#Low the hierarchy table that specifies the relationship between samples and predator species
sampleinfo <- read.csv("EuropeBatDiet.sampleinfo.csv",row.names=1)
#Visualise the whole table
sampleinfo
````

## Screening and editing the data files
We will first inspect the general properties of the OTU table.

````R
# Check how many samples are represented in the OTU table
ncol(otutable)

# Check how many OTUs are present in the OTU table
nrow(otutable)

# Check whether the OTU table has been normalised
colSums(otutable)

# We might also want to check what is the minimum representation of an OTU in a sample
min(otutable[otutable > 0])
# [1] 0.0001943025
# There is at least one OTU in a sample with a relative representation of 0.01% of the total reads.

# An histogram will give us a better overview of the distribution of the relative abundances. Note that 0 values have been skipped.
hist(otutable[otutable > 0])
````

At this point we might be interested in filtering OTUs with very low representation. We can do so using the copy_filt() function. This function can be used either for filtering OTUs based on raw abundance values (if using integers) or relative abundance values (if using decimals).
````R
# Check how many samples are represented in the OTU table
nrow(otutable)
otutable <- copy_filt(otutable,threshold=0.02)
````









Now, we will check whether the OTU names in the table and the tree match, as this is essential for the correct analysis of the data. This can be easily done using the function match_data().

````R
match_data(otutable,tree)
# The OTU tree contains OTUs absent in the OTU table. Filter the OTU tree.
````

The function warns that the OTU tree contains OTUs absent in the OTU table, and tells us to filter the OTU tree. The same function can be used to perform such a filtering, by specifying the type of output we want.
````R
match_data(otutable,tree,output="tree")
# The following OTUs/ASVs were removed from the tree for being absent in the count table: OTU334, OTU530, OTU1948, OTU1309, # OTU1970
# 
# Phylogenetic tree with 699 tips and 698 internal nodes.
# 
# Tip labels:
# 	OTU944, OTU2775, OTU168, OTU450, OTU700, OTU219, ...
# 
# Rooted; includes branch lengths.
````

The function returns that a few OTUs have been removed from the tree, but the tree has not been saved. To save the tree, we need to specify a new object, which can have the same name as the original tree. We can run the match_data() function again to ensure the OTU names at the OTU table and tree are matching.
````R
tree <- match_data(otutable,tree,output="tree")
match_data(otutable,tree)
# OTUs in the OTU table and OTU tree match perfectly.
````
