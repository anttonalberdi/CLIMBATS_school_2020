# Introduction to Hill numbers
Richness, Shannon index and Simpson index belong to a single statistical framework, as they all are monotonic functions of the basic sum, that is, the sum of the relative abundances of the types (pi) elevated to the q value (Jost, 2006; Keylock, 2005). This implies that Hill numbers (qD), or actual diversities, rather than entropies (e.g. Shannon index) or probabilities (e.g. Simpson index), can be formulated in terms of the basic sum (qλ) and the parameter q.

![Hill numbers expression](https://github.com/anttonalberdi/CLIMBATS_school_2020/blob/master/diversity_hillnumbers/images/hill-basic_sum.png)

This expression was first discovered by Hill (1973), hence the use of the name “Hill numbers” to refer to the output of this formula. Hill numbers have two major advantages over diversity indi‐ ces: (a) the interpretation of the measure and its measurement unit is always the same (Chao, Chiu, et al., 2014a; Tuomisto, 2010a), and ii) the sensitivity towards abundant and rare OTUs can be modulated with the parameter q. 

#### Effective number of OTUs
The expression yields a diver‐ sity measure in “effective number of OTUs”, that is, the number of equally abundant OTUs that would be needed to give the same value of diversity (Hill, 1973; Jost, 2006). When all OTUs in a system have the same relative abundances, as in the moth example given above, the effective number of OTUs for all q values equals the actual number of OTUs, namely richness.

````R
library(hilldiv)

# Create an even system composed of 30 OTUs each represented by 1 sequence
evensystem <- c(rep(1,30))
evensystem
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

# When the types are evenly distributed in the system, Hill numbers 
# of q-value 0, 1 and 2 (or any other q value) equal richness

index_div(evensystem,index="richness")
# [1] 30
hill_div(evensystem,qvalue=0)
# [1] 30
hill_div(evensystem,qvalue=1)
# [1] 30
hill_div(evensystem,qvalue=2)
# [1] 30
````

When the relative abundances of the types vary however, then the effective number of OTUs for q > 0 values decreases progresivelly. The higher the heterogeneity between types, the sharper will be decrease of the effective number of OTUs.

````R
library(hilldiv)

# Create an uneven system composed of 15 OTUs represented 
# by 1 sequence and 15 OTUs represented by 5 sequences
unevensystem <- c(rep(1,15),rep(5,15))
unevensystem
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5

# When types (OTUs) are not evenly distributed in the system, Hill numbers 
# of q-value > 0 show a decreasing effective number of types (OTUs)

index_div(unevensystem,index="richness")
# [1] 30
hill_div(unevensystem,qvalue=0)
# [1] 30
hill_div(unevensystem,qvalue=1)
# [1] 23.53793
hill_div(unevensystem,qvalue=2)
# [1] 20.76923
hill_div(unevensystem,qvalue=5)
# [1] 18.83793
````
In extreme cases in which the system is dominated by a few equally abundant OTUs, the effective number of OTUs will approach the number of those abundant OTUs.

````R
library(hilldiv)

# Create a super uneven system composed of 29 OTUs represented 
# by 1 sequence and 1 OTU represented by 971 sequences
superunevensystem <- c(rep(1,29),971)
superunevensystem
# [1]   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1
# [25]   1   1   1   1   1 971

# When types (OTUs) are very unevenly distributed in the system, Hill numbers 
# of q-value > 0 show a sharp drop of effective number of types (OTUs),
# approaching the number of dominant types (OTUs)

index_div(superunevensystem,index="richness")
# [1] 30
hill_div(superunevensystem,qvalue=0)
# [1] 30
hill_div(superunevensystem,qvalue=1)
# [1] 1.257225
hill_div(superunevensystem,qvalue=2)
# [1] 1.060592
hill_div(superunevensystem,qvalue=5)
# [1] 1.037471
````
#### Plotting diversity profiles
Hill numbers also enable diversity profiles of systems and subsystems to be plotted as continuous functions of the parameter q. This is useful to characterize the OTU abundance distri‐ bution of a system, as different compositions and abundance distributions can yield the same value for a particular order of diversity (e.g., q=1),but not for many of them(e.g. q=0, q=0.5 and q=1). Hill numbers convey all information contained in a species abundance distribution at a glance (Chao, Chiu, et al., 2014a; Leinster & Cobbold, 2012).

````R
library(hilldiv)

#Define the three model systems 
evensystem <- c(rep(1,30))
unevensystem <- c(rep(1,15),rep(5,15))
superunevensystem <- c(rep(1,29),971)

#Merge them in a single OTU table

mergedsystem <- cbind(evensystem,unevensystem,superunevensystem)
mergedsystem_profile <- div_profile(mergedsystem)
div_profile_plot(mergedsystem_profile)
````

#### Modulating the q value to adjust to study systems' properties
As shown above, the sensitivity towards abundant and rare OTUs can be modulated using the scaling parameter q, known as the “order” of diversity (Jost, 2006). **The larger the q value, the higher the importance attributed to abundant OTUs.** The ability to modulate the sensitivity towards abundant and rare OTUs by modifying a single parameter (q) is a useful means with which to adjust diversity measurements to the type of data and re‐ search question. For example, when rare types are considered to be of low importance (e.g., when attempting to define a core diet or mi‐ crobiome), or when rare types are considered untrustworthy due to technical issues (e.g., PCR or sequencing errors), researchers might opt for using a high q value, for example, q = 2, which overweighs abundant OTUs. The result can be interpreted as the effective num‐ ber of dominant OTUs in the system (Chao, Chiu, et al., 2014a). In contrast, if rare types are considered essential for the system, or re‐ searchers do not trust the relative abundance data due to potential technical biases, researchers might opt for using a q value of 0 that simply counts the number of types.

#### Relationship between Hill numbers and traditional popular indices
Three q values are particularly relevant, both for their significance, and their close relationship to popular diversity indices: q = 0, q = 1 and q = 2. Indeed, common diversity indices can be transformed to Hill numbers by applying simple mathematical transformations, as shown below.

As we have observed in the previous example, when a diversity of order zero (q = 0) is applied to the Hill numbers expression, it becomes insensitive to OTU frequencies, thus yielding a richness value. As the relative abundances of OTUs are overlooked, rare OTUs are overweighed.

````R
index_div(superunevensystem,index="richness")
# [1] 30
hill_div(superunevensystem,qvalue=0)
# [1] 30
````

A q value of 1 (in practical terms its limit, as the Hill number is undefined for q = 1) is the value that weighs OTUs by their frequency, without disproportionately favouring either rare or abundant ones (Jost, 2006). The value it yields is exactly the exponential of the Shannon index. In fact, q values under unity favour rare OTUs, while values above one favour abundant OTUs (Keylock, 2005). 

````R
hill_div(superunevensystem,qvalue=1)
# [1] 1.257225
index_div(superunevensystem,index="shannon")
# [1] 0.2289003
exp(index_div(superunevensystem,index="shannon"))
# [1] 1.257217
````

When a q value of 2 is applied, abundant OTUs are overweighed, and the formula yields the multiplicative inverse of the Simpson index.

````R
hill_div(superunevensystem,qvalue=2)
# [1] 1.060592
index_div(superunevensystem,index="simpson")
# [1] 0.05713
1/(1-index_div(superunevensystem,index="simpson"))
# [1] 1.060592
````
