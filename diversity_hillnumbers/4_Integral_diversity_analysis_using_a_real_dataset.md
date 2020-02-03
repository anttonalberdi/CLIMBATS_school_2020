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
