# Diversity partitioning and dissimilarity measurement
Dissimilarity and similarity metrics are based on partitioning of the diversity into different hierarchical values.

````R
library(hilldiv)

# Create an abundance table, with columns defining samples, and rows OTUs. 
abundance.table <- cbind(Sample1=c(0,0,10,10),Sample2=c(1,10,10,10),Sample3=c(10,10,0,0),Sample4=c(1,10,10,10))
rownames(abundance.table) <- c("OTU1","OTU2","OTU3","OTU4")
abundance.table 
#     Sample1 Sample2 Sample3 Sample4
# OTU1       0       1      10       1
# OTU2       0      10      10      10
# OTU3      10      10       0      10
# OTU4      10      10       0      10

````

````R
library(hilldiv)

#First, compute Hill numbers for each sample individually
hill_div(abundance.table,qvalue=0)
hill_div(abundance.table,qvalue=1)
hill_div(abundance.table,qvalue=2)
````

## Diversity partitioning
Diversity can also be computed for the whole system, and when doing so, there are different approaches that can be taken. In ecology, the idea of diversity has been traditionally broken down into three components: alpha (α), beta (β) and gamma (γ) diversities (Whittaker, 1960). In general terms, α‐diversity refers to the average diversity of subsystems or samples (although see below), β‐diversity measures the differences between subsystems (although see discussion below), while γ‐diversity includes the entire diversity of the system. Despite the existence of different approaches for diversity partitioning, within the framework of Hill number diversity partitioning responds to a multiplicative definition qDγ = qDα × qDß (Chao et al., 2012; Jost, 2007); that is, beta diversity is obtained by dividing gamma diversity by alpha diversity.

````R
#Compute alpha diversity values under different q values
alpha_div(abundance.table,qvalue=0)
alpha_div(abundance.table,qvalue=1)
alpha_div(abundance.table,qvalue=2)

#Compute gamma diversity values under different q values
gamma_div(abundance.table,qvalue=0)
gamma_div(abundance.table,qvalue=1)
gamma_div(abundance.table,qvalue=2)

#Compute beta diversity values under different q values
gamma_div(abundance.table,qvalue=0)/alpha_div(abundance.table,qvalue=0)
gamma_div(abundance.table,qvalue=1)/alpha_div(abundance.table,qvalue=1)
gamma_div(abundance.table,qvalue=2)/alpha_div(abundance.table,qvalue=2)
````

### Alpha diversity
It is important to note that alpha diversity is not obtained by averaging the Hill numbers of the subsystems, but computing the Hill numbers from the averaged basic sums of the subsystems (Chao et al., 2012). 

````R
# Alpha diversity
alpha_div(abundance.table,qvalue=1)
# [1] 2.584195

# Average of individual sample Hill numbers
mean(hill_div(abundance.table,qvalue=1))
# [1] 2.669514

# Averaged basic sums of the samples taken to the exponential (see relationship between 
# Shannon index and Hill number of q=1)
exp(mean(index_div(abundance.table,index="shannon")))
# [1] 2.584195
````

#### Unweighted
````R
alpha_div(abundance.table,qvalue=2)
# [1] 2.459373

#With equal weights
pi <- tss(abundance.table)
pi.q <- pi^qvalue
(sum(colSums(pi.q))/N)^(1/(1 - qvalue))
# [1] 2.459373
````

#### Weighted
````R
#With equal weights
weightvector <- rep(1/4,4)
qvalue=2

alpha_div(abundance.table,qvalue=2,weight=weightvector)
# [1] 2.459373

weight <- rep(1/4,4)
pi <- tss(abundance.table)
pi.w <- sweep(pi, 2, weight, "*")
pi.w.q <- pi.w^qvalue
N = ncol(abundance.table)
sum(rowSums(pi.w.q))^(1/(1 - qvalue))/N
# [1] 2.459373

#With unequal weights
weightvector <- c(0.2,0.2,0.3,0.3)
qvalue=2

alpha_div(abundance.table,qvalue=2,weight=weightvector)
# [1] 2.364782

pi <- tss(abundance.table)
pi.w <- sweep(pi, 2, weight, "*")
pi.w.q <- pi.w^qvalue
N = ncol(abundance.table)
sum(rowSums(pi.w.q))^(1/(1 - qvalue))/N
# [1] 2.364782
````

### Gamma diversity


````R
# Gamma diversity
gamma_div(abundance.table,qvalue=2)
# [1] 3.762173

weight <- rep(1/4,4)
pi <- tss(abundance.table)
pi.w <- sweep(pi, 2, weight, "*")
sum(rowSums(pi.w)^qvalue)^(1/(1 - qvalue))
# [1] 3.762173
````

### Beta diversity

Beta diversity is often used to vaguely refer to any kind of compositional heterogeneity among systems (Barwell, Isaac, & Kunin, 2015; Chao, Chiu, et al., 2014a; Tuomisto, 2010a,2010b). However, when diversity partitioning is carried out using Hill numbers, beta diversity is an actual diversity value that measures the **effective number of equally large and completely distinct subsystems** in a system. The Hill number beta diversity can also be interpreted as a unitless scalar that quantifies the **ratio of diversities between two levels (alpha and gamma)** of observation; thus, it also quantifies how many times richer an entire system is in effective OTUs (gamma diversity) than its constituent subsystems are on average (alpha diversity).


````R
#Diversity partitioning between Sample2 and Sample4 (identical systems)
div_part(abundance.table[,c(2,4)],qvalue=0)
# $Beta
# [1] 1

#Diversity partitioning between Sample1 and Sample3 (completely different systems)
div_part(abundance.table[,c(1,3)],qvalue=0)
# $Beta
# [1] 2

#Diversity partitioning between Sample1 and Sample2 (¡different systems)
div_part(abundance.table[,c(1,2)],qvalue=0)
# $Beta
# [1] 1.333333
````

## (Dis)similarity computation

````R
#Create a hierarchy table
hierarchy.table <- cbind(Samples=c("Sample1","Sample2","Sample3","Sample4"),Groups=c("SP1","SP1","SP2","SP2"))
hierarchy.table
#      Samples   Groups
# [1,] "Sample1" "SP1" 
# [2,] "Sample2" "SP1" 
# [3,] "Sample3" "SP2" 
# [4,] "Sample4" "SP2"

#Diversity partitioning between Sample2 and Sample4
div_part(abundance.table,hierarchy=hierarchy.tableπ[c(2,4),])


````
