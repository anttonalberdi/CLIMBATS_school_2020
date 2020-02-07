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
Diversity can also be computed for the whole system, and when doing so, there are different approaches that can be taken. In ecology, the idea of diversity has been traditionally broken down into three components: alpha (α), beta (β) and gamma (γ) diversities (Whittaker, 1960). In general terms, α‐diversity refers to the average diversity of subsystems or samples (although see discussion about the different α‐diversities in Chao, Chiu, & Hsieh, 2012), β‐diversity measures the differences between subsys‐ tems (although see discussion about the different β‐diversities in Tuomisto, 2010a), while γ‐diversity includes the entire diversity of the system (Figure 5). Despite the existence of different ap‐ proaches for diversity partitioning, within the framework of Hill number diversity partitioning responds to a multiplicative defini‐ tion qDγ = qDα × qDß (Chao et al., 2012; Jost, 2007); that is, beta diversity is obtained by dividing gamma diversity by alpha diversity.

````R
#Compute alpha diversity values under different q values
alpha_div(abundance.table,qvalue=0)
alpha_div(abundance.table,qvalue=1)
alpha_div(abundance.table,qvalue=2)

#Compute gamma diversity values under different q values
gamma_div(abundance.table,qvalue=0)
gamma_div(abundance.table,qvalue=1)
gamma_div(abundance.table,qvalue=2)
````

### Alpha diversity

````R

````

### Gamma diversity

### Beta diversity

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
````