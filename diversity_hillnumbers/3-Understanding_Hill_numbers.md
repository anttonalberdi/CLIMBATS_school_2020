# 3 - Understanding Hill numbers 
The aim of this lesson is to get familiar with the Hill numbers expression and understand how and why the diversity values vary when using systems with different relative representation of OTUs, different q-values and differently-shaped trees. 

![Hill numbers expression](https://github.com/anttonalberdi/CLIMBATS_school_2020/blob/master/diversity_hillnumbers/images/hill-basic_sum.png)

For doing so, and for the sake of simplicity, we will work with mock systems rather than OTU tables. The systems we will create will be comprised of 30 OTUs.

````R
library(hilldiv)

evensystem <- c(rep(1,30))
names(evensystem) <- paste("OTU",c(1:30),sep="")

unevensystem <- c(rep(1,15),rep(5,15))
names(unevensystem) <- paste("OTU",c(1:30),sep="")

superunevensystem <- c(rep(1,29),971)
names(superunevensystem) <- paste("OTU",c(1:30),sep="")
````

## Neutral Hill numbers (without considering phylogenies)

### Even system

We will start with the simplest of the systems. The diversity of the system for all q values is the same, namely the number of OTUs, because the relative representation of each OTU is the identical.

````R
hill_div(evensystem,qvalue=0)
# [1] 30
hill_div(evensystem,qvalue=1)
# [1] 30
hill_div(evensystem,qvalue=2)
# [1] 30
````

The mathematical procedure to reach that diversity value is different though.

````R
vector <- evensystem
pi <- vector[vector!=0]
# pi
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
sum(pi)
# 30
````

The Hill numbers expression requires each OTU to be represented in relative rathern than absolute terms, i.e. the sum of all OTU representation values must sum to 1. Be aware that all hilldiv functions include this step normalisation step. In this example we need to do it ourselves. The hilldiv function tss() will help us to normalise the data very easily.

````R
pi <- tss(vector[vector!=0])
# [1] 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333
# [14] 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333
# [27] 0.03333333 0.03333333 0.03333333 0.03333333
sum(pi)
# [1] 1
````

Next, we will apply the Hill numbers expression using different q values.

````R
qvalue=0
sum(pi^qvalue)^(1/(1-qvalue))
# [1] 30
````

As expected, the diversity value is identical to what "hill_div(evensystem,qvalue=0)" yielded:

````R
qvalue=1
sum(pi^qvalue)^(1/(1-qvalue))
# [1] 1
````

When using the q-value 1 the result is not identical to what we got with "hill_div(evensystem,qvalue=1)". The reason is that the Hill number is undefined for q=1. We won't explain why in this lesson. The limits of unity are defined though, which yield an almost identical results: 

````R
qvalue=0.99999999
sum(pi^qvalue)^(1/(1-qvalue))
# [1] 30
````

And finally the same for q=2

````R
qvalue=2
sum(pi^qvalue)^(1/(1-qvalue))
# [1] 30
````
Next, we will split the expression "sum(pi^qvalue)^(1/(1-qvalue))" in two, to observe how the q-value affects the relative representation of the OTUs, the exponent and the final result.

````R
# When q=0
qvalue=0
pi^qvalue
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

1/(1-qvalue)
# 1

sum(pi^qvalue)^(1/(1-qvalue))
# 30

# When q=2
qvalue=0
pi^qvalue
# [1] 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111
# [13] 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111
# [25] 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111

1/(1-qvalue)
[1] -1

sum(pi^qvalue)^(1/(1-qvalue))
# 30


````

````R
hill_div(superunevensystem,qvalue=0)
# [1] 30
hill_div(superunevensystem,qvalue=1)
# [1] 1.257225
hill_div(superunevensystem,qvalue=2)
# [1] 1.060592
````

Next, we can do the same operations with the super uneven system

````R
vector <- superunevensystem
pi <- tss(vector[vector!=0])
pi
# [1] 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
# [26] 0.001 0.001 0.001 0.001 0.971

qvalue=0
pi^qvalue
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

1/(1-qvalue)
# [1] 1

sum(pi^qvalue)^(1/(1-qvalue))
# [1] 30

qvalue=2
pi^qvalue
# [1] 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001
# [17] 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.942841

1/(1-qvalue)
# [1] -1

sum(pi^qvalue)^(1/(1-qvalue))
# [1] 1.060592
````

## Phylogenetic Hill numbers

The mathematical operations become slightly more complicated when including trees into the computation of diversities. The effect of different-shaped trees can be clearly shown though when comparing a simple star-shaped tree with a hierarchical tree. First of all, let's create those simple trees and check how the Hill number value is changed.

````R
#Create an even system
evensystem <- c(5,5,5)
names(evensystem) <- c("OTU1","OTU2","OTU3")

#Create trees
library(phytools)
eventree <- starTree(c("OTU1","OTU2","OTU3"), branch.lengths=c(1,1,1))
uneventree <- read.tree(text="((OTU1:0.5,OTU2:0.5):0.5,OTU3:1);")

hill_div(evensystem,qvalue=0)
[1] 3
hill_div(evensystem,qvalue=0,tree=eventree)
[1] 3
hill_div(evensystem,qvalue=0,tree=uneventree)
[1] 2.5
````
![Example trees](https://github.com/anttonalberdi/CLIMBATS_school_2020/blob/master/diversity_hillnumbers/images/example_trees.png)

The even tree yields exactly the same diversity value as the neutral Hill number, while the hierarchical tree yields a lower diversity value. Remember that while the neutral Hill number yields a diversity value in "effective number of OTUs", the phylogenetic Hill number yields a diversity value in "effective number of lineages".

The computation of phylogenetic Hill numbers requires some operations to be done before computing diversity:
````R
#All these are preparations
Li <- tree$edge.length
ltips <- sapply(tree$edge[, 2], function(node) tips(tree, node))
ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
T <- sum(Li * ai)
Li <- Li[ai != 0]
ai <- ai[ai != 0]

#To apply the phylogenetic Hill number formula
sum(Li/T * ai^qvalue)^(1/(1-qvalue))
````
In the following, we will analyse each step one by one.

### Even tree

````R
vector <- tss(c(5,5,5))
names(vector) <- c("OTU1","OTU2","OTU3")
tree <- eventree
qvalue=0

Li <- tree$edge.length
Li
# [1] 1 1 1

library(geiger)
ltips <- sapply(tree$edge[, 2], function(node) tips(tree, node))
ltips
# [1] "OTU1" "OTU2" "OTU3"

ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
ai
# [1] 0.3333333 0.3333333 0.3333333

T <- sum(Li * ai)
T
# [1] 1

Li <- Li[ai != 0]
Li
# [1] 1 1 1

ai <- ai[ai != 0]
ai
# [1] 0.3333333 0.3333333 0.3333333

sum(Li/T * ai^qvalue)^(1/(1-qvalue))
# 3

Li/T
# [1] 1 1 1

ai^qvalue
# [1] 0.3333333 0.3333333 0.3333333

1/(1-qvalue)
# [1] 1

````

### Uneven tree

````R
vector <- tss(c(5,5,5))
names(vector) <- c("OTU1","OTU2","OTU3")
tree <- uneventree
qvalue=0

vector
#      OTU1      OTU2      OTU3 
# 0.3333333 0.3333333 0.3333333 

Li <- tree$edge.length
Li
# [1] 0.5 0.5 0.5 1.0

library(geiger)
ltips <- sapply(tree$edge[, 2], function(node) tips(tree, node))
ltips
# [[1]]
# [1] "OTU1" "OTU2"

# [[2]]
# [1] "OTU1"

# [[3]]
# [1] "OTU2"

# [[4]]
# [1] "OTU3"

ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
ai
# [1] 0.6666667 0.3333333 0.3333333 0.3333333

T <- sum(Li * ai)
T
# [1] 1

Li <- Li[ai != 0]
Li
# [1] 0.5 0.5 0.5 1.0

ai <- ai[ai != 0]
ai
# [1] 0.6666667 0.3333333 0.3333333 0.3333333

sum(Li/T * ai^qvalue)^(1/(1-qvalue))
# 2.5

Li/T
# [1] 0.5 0.5 0.5 1.0

ai^qvalue
# [1] 1 1 1 1

1/(1-qvalue)
# [1] 1

````
