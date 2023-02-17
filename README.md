# Resolving Higher Level Networks using a polytypic salamander system

Code from the paper “Resolving higher-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Plethodontidae: Desmognathus)”

### Overview
Code in this repo is from the paper “Resolving higher-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Plethodontidae: Desmognathus)”.

We include two examples pipelines, one one simulated data that demonstrates a proof of concept, and the other used on our empirical data. We take each in turn below. 

### Simulated Data
The directory called `SimulatedData_analysis` has a single R script called `Obtain_DisplayedTrees_Example.R` which has the steps to obtain the displayed trees for the data simulated via Hybrid-Lambda on the network:
```
	 "(((1:2,h_2#.4:1.5)R:3,((((2:.5)h_2#.4:1,4:1.5)ff:.5,(3:1,(5:.5)h_1#.7:.5)g:1)e:1,((h_1#.7:.5, 6:1)h:1,7:2)f:1)c:2)a:2)r;"
```
We divide our method into several steps. Steps 1 and 2 generate a quartet table using Nanuq. Step 3 takes the output from Step 2 and classifies all quartets as either non-deep or deep-hybrid and shows the tree (or trees) associated with such quartets. Step 4 constructs all trees that display trees associated with non-deep quartet networks, and one of the quartet trees associated to each quartet identified as deep-hybrid. Step 5 identifies missing quartets displayed in the network for each tree found in Step 4, and step 6 tallies each type of quartet and displays the estimated displayed trees. Steps 3–6 rely on functions within `DeepQuartets_functions.R`:
- Quartet_Analysis: This function performs Step (3) of our method. 
- Trees_quartet_display: This function performs Step (4) of our method. 
- Missing_quartets_on_trees: This function performs Step (5) of our method. 
- Tallying_quartets_on_trees: This function performs Step (6) of our method. 

This method assumes you have 7 taxa, and as such comes with an R data file of all possible 7-taxa trees `7taxaTrees.RData` (needed to compare all 7-taxon topologies exhaustively) as well as simulated gene trees for 7 taxa `Example_network_sample`.

Finally, make sure before you 'source' the functions R script, you open that file in a text editor and update the wd path.

### Anchored Hybrid Enrichment Data
+ The AHE pipeline starts with `1_Sample_trees_AHE.R` which contains the routine for subsampling gene tree from the AHE data set using the functions in the file 'DeepQuartets_functions.R'. Specifically, we select one individual at random for each species. Then we keep a copy of each gene tree that displays these individuals and removed all other individuals. This left a set of gene trees with 7 taxa each. This is repeated 10,000 times.
+ `Obtain_DisplayedTrees_AHE.R`: This file contains the steps to obtain the backbone tree for the AHE dataset. This whole method was outlined in the section above. 
+ `UseDeepHybrid_SNF_Dsuite.R`: This file contains the integration of the information of sNMF and Dsuite and the Deep quartet network information. Here we are just visualizing the integration of admixture signals across our network and other analyses from our study to create a final higher-level network.


