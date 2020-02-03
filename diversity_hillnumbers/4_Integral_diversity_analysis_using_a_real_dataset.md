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

## Screening the data files
We will first inspect the general properties of the OTU table and OTU tree.

````R
#Check how many samples are represented in the OTU table
ncol(otutable)

#Check how many OTUs are present in the OTU table
nrow(otutable)

#Check whether the OTU table has been normalised
colSums(otutable)

#Obtain a general overview of the tree properties
tree
````

Now, we will check whether the OTU names in the table and the tree match, as this is essential for the correct analysis of the data. This can be easily done using the function match_data().

````R
match_data(otutable,tree)
# The OTU tree contains OTUs absent in the OTU table. Filter the OTU tree.
````

The function warns that the OTU tree contains OTUs absent in the OTU table, and tells us to filter the OTU tree. The same function can be used to perform such a filtering, by specifying the type of output we want.
````R
match_data(otutable,tree,output="tree")
````

The function returns that a few OTUs have been removed from the tree, but the tree has not been saved. To save the tree, we need to specify a new object, which can have the same name as the original tree. We can run the match_data() function again to ensure the OTU names at the OTU table and tree are matching.
````R
tree <- match_data(otutable,tree,output="tree")
match_data(otutable,tree)
# OTUs in the OTU table and OTU tree match perfectly.
````
