# Resolving Higher Level Networks using a polytypic salamander system

Code from the paper “Resolving higher-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Plethodontidae: Desmognathus)”

### Overview
The code in this repo is from the paper “Resolving higher-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Plethodontidae: Desmognathus)”. We briefly give an overview of the steps of the method and their implementation. All the functions mentioned in these steps can be found in the file DeepQuartets_functions.R unless it is stated otherwise.

##### Step 1
The first step is determining the empirical quartet counts (empirical concordance factors) from the gene trees for every subset of 4 taxa. This step was previously implemented in the R package 'MSCquartets' in the function NANUQ() by Allman et al. (2020) (whose input and output are discussed in the next step). We used such an implementation to obtain all empirical quartet counts. 

##### Step 2
The second step is to apply a statistical hypothesis test to the empirical concordance factors (obtained from Step 1) as done by Mitchell et al. (2019) to determine which quartet networks are B-quartet networks and which are T-quartet networks.
The null hypothesis is that the quartet is a T-quartet, thus if we reject the null hypothesis we consider the quartet as a B-quartet. This has also been implemented in the R package 'MSCquartets' in the function NANUQ(). The input of such function input is a collection of gene trees and the output is a table with the empirical concordance factors with the respective p-values. We used such implementation to determine which quartets are B-quartets and which are T-quartets.

##### Step 3
The third step is to associate one quartet tree with each quartet identified as a T-quartet and associate two quartet trees with each quartet identified as a B-quartet. As explained in the manuscript, we expect that the collection of all quartet trees associated with B-quartets and T-quartets will be displayed in the network. We refer to this collection as the set of dominant quartet trees. As detailed in the SI, for each T-quartet the dominant quartet tree is the most frequent one and for each B-quartet the dominant quartets are the two most frequent quartet trees. 

For this step, we wrote the R function Quartet_Analysis(). The input for this function is the table outputted by NANUQ() and a significance level to determine whether we accept or reject the null hypothesis (as explained in the previous step). The output is a list with the dominant quartets obtained from T-quartets (1 form each), and a list with the dominant quartets obtained from all B-quartets (two from each). 

##### Step 4
The fourth step consists of constructing a set T^0 consisting of all trees that display only dominant quartet trees (those quartet trees in step 3). Any tree in T0 could potentially be displayed in the network since, for any subset of 4 taxa, each of these trees will display a dominant quartet tree on those 4 taxa. For now, this step is done exhaustively by looking at all possible trees displaying only dominant quartets. This is possible since for our data set and simulated data sets, we only consider 7 taxa. To obtain such a set of trees we implemented a function called Trees_quartet_display() whose input is the list of dominant quartets as output by the function Quartet_Analysis() explained in the previous step, and a vector with the name of the taxa of the seven taxa in the data. This function finds all trees on 7 taxa whose quartet trees are all a subset of the set of dominant quartets. 

##### Step 5
The fifth step consists of, for each tree T in T^0, constructing a set Q(T) consisting of the set of quartet trees associated with all B-quartets not displayed in T. Note that for any candidate tree T in T^0, Q(T) contains all the dominant quartet trees that are not displayed in T. If one would assume that T is a displayed tree in the network, then there should be other displayed trees in the network that display the remaining dominant quartet trees, i.e. those in Q(T). Loosely speaking, Q(T) is the information in the network not accounted for by T.  For this step we implemented the function Missing_quartets_on_trees(). The input is the set of trees as produced by the function Trees_quartet_display() explained in the previous step. This function compares each tree with the set of dominant trees and for each tree T, it outputs a list of the dominant quartets not displayed in T (which is referred to in the paper as Q(T)).

##### Step 6
For each tree T in T^0, find the set of trees S(T) contained in T^0 that display the maximum number of quartet trees in Q(T). The output of the method consists of those trees in T^0 such that the trees in S(T) displayed the greatest number of quartets for all T in T^0, and their associated trees are the set S(T). We refer to S(T) as the complement of T. In particular, note that for any candidate tree T in T^0, the set S(T) contains the trees that display most of the dominant quartet trees not displayed by T. One can think of S(T) as the set of candidate trees that display most of the information missing in T. The output consists of those trees T whose complement describes most of the missing data in T, together with S(T). For this step, we implemented the function Tallying_quartets_on_trees() whose input is the set of trees T^0 as output by Trees_quartet_display(), the list of quartets not displayed by a tree (Q(T)), as output by the function Missing_quartets_on_trees() and the set of taxon names. The output consists of a collection of trees displayed in the network.

##### Data
We include two example pipelines, one one simulated data that demonstrates a proof of concept, and the other used on our empirical data. We take each in turn below. 

### Simulated Data
The directory called `SimulatedData_analysis` has a single R script called `Obtain_DisplayedTrees_Example.R` which has the steps to obtain the displayed trees for the data simulated via Hybrid-Lambda on the network:
```
	 "(((1:2,h_2#.4:1.5)R:3,((((2:.5)h_2#.4:1,4:1.5)ff:.5,(3:1,(5:.5)h_1#.7:.5)g:1)e:1,((h_1#.7:.5, 6:1)h:1,7:2)f:1)c:2)a:2)r;"
```
We divide our method into several steps. Steps 1 and 2 generate a quartet table using Nanuq. Step 3 takes the output from Step 2, classifies all quartets as non-deep or deep-hybrid, and shows the tree (or trees) associated with such quartets. Step 4 constructs all trees that display trees associated with non-deep quartet networks, and one of the quartet trees associated to each quartet identified as deep-hybrid. For each tree found in Step 4, Step 5 identifies quartets displayed in the network not displayed in the tree, and Step 6 finds the trees that display those quartets found in Step 5 and displays the estimated displayed trees. Steps 3–6 rely on functions within `DeepQuartets_functions.R`:
- Quartet_Analysis: This function performs Step (3) of our method. 
- Trees_quartet_display: This function performs Step (4) of our method. 
- Missing_quartets_on_trees: This function performs Step (5) of our method. 
- Tallying_quartets_on_trees: This function performs Step (6) of our method. 

This method assumes you have 7 taxa, and as such comes with an R data file of all possible 7-taxa trees `7taxaTrees.RData`. This file is needed to be able to compare all 7-taxon topologies exhaustively. We also include simulated gene trees on 7 taxa `Example_network_sample`.

Finally, make sure before you 'source' the functions R script, you open that file in a text editor and update the wd path.

### Anchored Hybrid Enrichment Data
+ The AHE pipeline starts with `1_Sample_trees_AHE.R` which contains the routine for subsampling gene tree from the AHE data set using the functions in the file 'DeepQuartets_functions.R'. Specifically, we select one individual at random for each species. Then we keep a copy of each gene tree that displays these individuals and removed all other individuals. This left a set of gene trees with 7 taxa each. This is repeated 10,000 times.
+ `Obtain_DisplayedTrees_AHE.R`: This file contains the steps to obtain the backbone tree for the AHE dataset. This whole method was outlined in the section above. 
+ `UseDeepHybrid_SNF_Dsuite.R`: This file contains the integration of the information of sNMF and Dsuite and the Deep quartet network information. Here we are just visualizing the integration of admixture signals across our network and other analyses from our study to create a final higher-level network.


