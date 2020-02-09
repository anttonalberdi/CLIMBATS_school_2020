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

## Screening and pre-processing the data files
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
# [1] 0.0001916443
# There is at least one OTU in a sample with a relative representation of 0.1% of the total reads.

# An histogram will give us a better overview of the distribution of the relative abundances. Note that 0 values have been skipped.
hist(otutable[otutable > 0])
````

At this point we might be interested in filtering OTUs with very low representation. We can do so using the copy_filt() function. This function can be used either for filtering OTUs based on raw abundance values (if using integers) or relative abundance values (if using decimals).
````R
# Check how many samples are represented in the OTU table
nrow(otutable)
# [1] 843
otutable <- copy_filt(otutable,threshold=0.001)
nrow(otutable)
# [1] 439
````

Now that we have removed a large quantity of OTUs, we will check whether the OTU names in the table and the tree match, as this is essential for the correct analysis of the data. This can be easily done using the function match_data().

````R
match_data(otutable,tree)
# The OTU tree contains OTUs absent in the OTU table. Filter the OTU tree.
````

As expected, the function warns that the OTU tree contains OTUs absent in the OTU table, and tells us to filter the OTU tree. The same function can be used to perform such a filtering, by specifying the type of output we want.

````R
match_data(otutable,tree,output="tree")
# The following OTUs/ASVs were removed from the tree for being absent in the count table: OTU1948, OTU1914, OTU615, OTU2578 (...)
# 
# Phylogenetic tree with 439 tips and 438 internal nodes.
# 
# Tip labels:
#	OTU168, OTU1787, OTU786, OTU1509, OTU358, OTU219, ...
# 
# Rooted; includes branch lengths.
````

The function returns that a few OTUs have been removed from the tree, but the tree has not been saved. To save the tree, we need to specify a new object, which can have the same name as the original tree. We can run the match_data() function again to ensure the OTU names at the OTU table and tree are matching.

````R
tree <- match_data(otutable,tree,output="tree")
match_data(otutable,tree)
# OTUs in the OTU table and OTU tree match perfectly.
````
## Beginning with the diversity analyses
We will first compute the Hill numbers of all samples independently.

````R
hill_div(otutable,qvalue=0)
hill_div(otutable,qvalue=1)
hill_div(otutable,qvalue=2)
````
That was fast. Remember we can also compute traditional diversity indices:

````R
index_div(otutable,index="richness")
index_div(otutable,index="shannon")
index_div(otutable,index="simpson")
````
Now the same for the phylogenetic Hill numbers.

````R
hill_div(otutable,qvalue=0,tree=tree)
hill_div(otutable,qvalue=1,tree=tree)
hill_div(otutable,qvalue=2,tree=tree)
````
Computing phylogenetic Hill numbers takes considerably longer than computing the neutral diversity metrics. This is because entropy values are calculated for every branch in the phylogenetic tree. In this example, the tree only has 439 tips and 438 internal nodes, so the computation is fast. When trees contain thousands of OTUs though, computing entropy values for all branches might take hours.

## Comparing sample-level diversity means across species
One of the basic operations when studying diversity patterns is to compare diversities across groups. In our example, we have 40 samples belonging to 4 different predator species. The relation between both is set by the so-called hierarchy table, in this example identified as "sampleinfo".

````R
#Visualise the structure
sampleinfo

#Test whether the mean diversity values differ across species when q=0
div_test(otutable,qvalue=0,hierarchy=sampleinfo)

#Test whether the mean diversity values differ across species when q=1
div_test(otutable,qvalue=1,hierarchy=sampleinfo)

#The same considering phylogenetic relations across OTUS
div_test(otutable,qvalue=0,hierarchy=sampleinfo,tree=tree)
div_test(otutable,qvalue=1,hierarchy=sampleinfo,tree=tree)
````
Short explanation

````R
#Create div_test objects
divq0 <- div_test(otutable,qvalue=0,hierarchy=sampleinfo)
divq1 <- div_test(otutable,qvalue=1,hierarchy=sampleinfo)
divq0phylo <- div_test(otutable,qvalue=0,hierarchy=sampleinfo,tree=tree)
divq1phylo <- div_test(otutable,qvalue=1,hierarchy=sampleinfo,tree=tree)

#Plot them using div_test_plot()
div_test_plot(divq0)
div_test_plot(divq1phylo, chart="jitter")
````
Running tests with post-hoc analyses

````R
divq0 <- div_test(otutable,qvalue=0,hierarchy=sampleinfo,posthoc=TRUE)
divq1 <- div_test(otutable,qvalue=1,hierarchy=sampleinfo,posthoc=TRUE)
divq0phylo <- div_test(otutable,qvalue=0,hierarchy=sampleinfo,tree=tree,posthoc=TRUE)
divq1phylo <- div_test(otutable,qvalue=1,hierarchy=sampleinfo,tree=tree,posthoc=TRUE)
````

Plot posthoc information

````R
#All p-values
div_test_plot(divq1phylo, chart="jitter",posthoc=TRUE)

#Significant p-values
div_test_plot(divq1phylo, chart="jitter",posthoc=TRUE,threshold=0.05)
````

Safe plots as PDF

````R
pdf("divq0.pdf",width=8,height=6)
div_test_plot(divq0, chart="jitter",posthoc=TRUE,threshold=0.05)
dev.off()

pdf("divq1.pdf",width=8,height=6)
div_test_plot(divq1, chart="jitter",posthoc=TRUE,threshold=0.05)
dev.off()

pdf("divq1phylo.pdf",width=8,height=6)
div_test_plot(divq1phylo, chart="jitter",posthoc=TRUE,threshold=0.05)
dev.off()
````

Ecological explanation
