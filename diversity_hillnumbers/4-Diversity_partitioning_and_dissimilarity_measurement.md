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

The Hill number beta diversity always ranges from 1 (when all subsystems are identical) to the actual number of subsystems (when all subsystems are completely different) (Chao, Chiu, et al., 2014a; Chiu et al., 2014).

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
Dissimilarity indices range between 0 and 1; 0 indicates that the subsystems compared are identical, while 1 indicates that they are completely different. As the beta diversity lies in between 1 and the total number of subsystems, the Hill number beta diversity cannot directly be used to compute dissimilarities. However, it is possible —and desirable— to remove the dependence on the number of subsystems and compute dissimilarity measures by applying simple transformations to beta diversity, both for diversities (Chao et al., 2012; Jost, 2007) as well as phylodiversities (Chao, Chiu, et al., 2014a; Chiu et al., 2014). 

Four classes of **similarity** measures derived from Hill number beta diversities have been proposed, from which dissimilarity measures can be obtained by calculating their one‐complements (1‐XqN). The Sørensen‐type classes quantify similarity from the perspective of the subsystem, while the Jaccard‐type classes quantify similarity from the perspective of the overall system (Chao et al., 2019; Chiu et al., 2014).

### Sørensen‐type overlap (CqN)
The Sørensen‐type overlap (CqN) quantifies the effective average proportion of a sub‐ system's OTUs (or lineages in the case of phylodiversities) that is shared across all subsystems. This is thus a metric that quantifies overlap from the subsystem's perspective. Its corresponding dissimilarity measure (1 − CqN) quantifies the effective average proportion of nonshared OTUs or lineages in a system.

````R
beta <- div_part(abundance.table[,c(1,2)],qvalue=0)$Beta
CqN(beta,qvalue=0,N=2)
# [1] 0.6666667

beta <- div_part(abundance.table[,c(1,2)],qvalue=1)$Beta
CqN(beta,qvalue=1,N=2)
# [1] 0.7947585

beta <- div_part(abundance.table[,c(1,2)],qvalue=2)$Beta
CqN(beta,qvalue=2,N=2)
# [1] 0.7933461
````

### Jaccard‐type overlap (CqN)
The Jaccard‐type overlap (UqN) quantifies the effective proportion of OTUs or lineages in a system that are shared across all subsystems. Hence, this metric quantifies overlap from the perspective of the overall system. Its corresponding dissimilarity (1 − UqN) quantifies the effective proportion of nonshared OTUs or lineages in the overall system.

````R
beta <- div_part(abundance.table[,c(1,2)],qvalue=0)$Beta
UqN(beta,qvalue=0,N=2)
# [1] 0.5

beta <- div_part(abundance.table[,c(1,2)],qvalue=1)$Beta
UqN(beta,qvalue=1,N=2)
# [1] 0.7947574

beta <- div_part(abundance.table[,c(1,2)],qvalue=2)$Beta
UqN(beta,qvalue=2,N=2)
# [1] 0.8847663
````

### Sørensen‐type turnover-complement (VqN)
The Sørensen‐type turnover‐complement (Vqn) is the complement of the Sørensen‐type turnover, which quantifies the normalized OTU turnover rate with respect to the average subsystem (i.e., alpha), thus provides the proportion of a typical subsystem that changes across subsystems (Harrison, Ross, & Lawton, 1992; Jost, 2007).

````R
beta <- div_part(abundance.table[,c(1,2)],qvalue=0)$Beta
VqN(beta,N=2)
# [1] 0.6666667

beta <- div_part(abundance.table[,c(1,2)],qvalue=1)$Beta
VqN(beta,N=2)
# [1] 0.8471202

beta <- div_part(abundance.table[,c(1,2)],qvalue=2)$Beta
VqN(beta,N=2)
# [1] 0.8847663
````

### Jaccard‐type turnover-complement (SqN)
The Jaccard‐type turnover‐complement (SqN)is the complement of the Jaccard‐type turnover, which quantifies the normalized OTU turnover rate with respect to the whole system (i.e. gamma).

````R
beta <- div_part(abundance.table[,c(1,2)],qvalue=0)$Beta
SqN(beta,N=2)
# [1] 0.5

beta <- div_part(abundance.table[,c(1,2)],qvalue=1)$Beta
SqN(beta,N=2)
# [1] 0.7347863

beta <- div_part(abundance.table[,c(1,2)],qvalue=2)$Beta
SqN(beta,N=2)
# [1] 0.7933461
````

These (dis)similarity metrics are generalisations for multiple systems and q-values that encompass, as special cases, some of the most popular (dis)similarity measures used in ecology (Chao, Chiu, et al., 2014a; Chao et al., 2016; Jost, 2007). For instance, C02 (the Sørensen‐type overlap between two systems [N = 2] when OTU phylogenies are not considered and q = 0) produces the Sørensen similarity index. 

````R
beta <- div_part(abundance.table[,c(1,2)],qvalue=0)$Beta
CqN(beta,qvalue=0,N=2)
# [1] 0.6666667

library(vegan)
betadiver(x=t(abundance.table[,c(1,2)]), method = "sor")
#           Sample1
# Sample2 0.6666667
````
Another noteworthy example is that the measure 1 − U02 (the one‐complement of the Jaccard‐type overlap when OTU phylogenies are considered, q = 0 and N = 2) is identical to the unweighted UniFrac distance (Lozupone & Knight, 2005):

````R
uneventree <- read.tree(text="(((OTU1:0.5,OTU2:0.5):0.25,OTU4:0.75):0.25,OTU3:1);")

unifrac <- GUniFrac(t(abundance.table[,c(1,2)]), uneventree)$unifracs
unifrac[, , "d_UW"]
#          Sample1   Sample2
# Sample1 0.0000000 0.3846154
# Sample2 0.3846154 0.0000000
 
beta <- div_part(abundance.table[,c(1,2)],qvalue=0,tree=uneventree)$Beta
1-UqN(beta,qvalue=0,N=2)
# [1] 0.3846154
````
Further relations between these four (dis)similarity measures and other popular in‐ dices can be found elsewhere (e.g., Jost, 2007, Chao et al., 2012, Chiu et al., 2014). If researchers opt for basing diversity measurements on Hill numbers, it is also advisable to frame dissimilarity measurements within the same scheme. Basing dissimilarity measurements on beta diversities derived from Hill numbers enables logical consistency to be kept with the conclusions based on Hill numbers (Chao et al., 2012; Jost, 2007). Furthermore, as all measures are continuous as q ranges from zero to infinity, (dis)similarity profiles can be made for any of them (Chiu et al., 2014).
