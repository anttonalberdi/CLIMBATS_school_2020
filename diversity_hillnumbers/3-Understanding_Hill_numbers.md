# 3 - Understanding Hill numbers 

````R
library(hilldiv)

evensystem <- c(rep(1,30))
unevensystem <- c(rep(1,15),rep(5,15))
superunevensystem <- c(rep(1,29),971)
````

## Neutral Hill numbers (without considering phylogenies)

````R
hill_div(evensystem,qvalue=0)
# [1] 30
hill_div(evensystem,qvalue=1)
# [1] 30
hill_div(evensystem,qvalue=2)
# [1] 30
````

````R
vector <- evensystem
pi <- vector[vector!=0]
# pi
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

pi <- tss(vector[vector!=0])
# [1] 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333
# [14] 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333
# [27] 0.03333333 0.03333333 0.03333333 0.03333333

qvalue=0
sum(pi^qvalue)^(1/(1-qvalue))
# [1] 30

qvalue=1
sum(pi^qvalue)^(1/(1-qvalue))
# [1] 1

qvalue=0.99999999
sum(pi^qvalue)^(1/(1-qvalue))
# [1] 30

qvalue=2
sum(pi^qvalue)^(1/(1-qvalue))
# [1] 30

pi^qvalue
# [1] 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111
# [13] 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111
# [25] 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111 0.001111111

1/(1-qvalue)
[1] -1

qvalue=0
pi^qvalue
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

1/(1-qvalue)
# 1
````

````R
hill_div(superunevensystem,qvalue=0)
# [1] 30
hill_div(superunevensystem,qvalue=1)
# [1] 1.257225
hill_div(superunevensystem,qvalue=2)
# [1] 1.060592
````

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
# 1

qvalue=2
pi^qvalue
# [1] 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001
# [17] 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.000001 0.942841

1/(1-qvalue)
# [1] -1
````

## Phylogenetic Hill numbers

````R
#Create an even system
evensystem <- c(5,5,5)
names(evensystem) <- c("OTU1","OTU2","OTU3")

#Create trees
library(phytools)
eventree <- starTree(c("OTU1","OTU2","OTU3"), branch.lengths=c(1,1,1))
uneventree <- read.tree(text="((OTU1:0.5,OTU2:0.5):0.5,OTU3:1);")
````

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
