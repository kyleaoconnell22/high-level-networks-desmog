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
We divide our method into several steps. Steps 1 and 2 are to generate a quartet table using Nanuq. Step 3 takes the output from Step 2 and classifies all quartets as wither non-deep or deep-hybrid and shows the tree (or trees) associated to such quartets. Step 4 constructs all trees that display trees associated with non-deep quartet networks, and one of the quartet trees associated to each quartet identified as deep-hybrid. Step 5 identifies missing quartets on the trees, and step 6 displays the estimated backbone tree. Steps 3–6 rely on functions within `DeepQuartets_functions.R`:
- Quartet_Analysis: This function performs Step (3) of our method. 
- Trees_quartet_display: This function performs Step (4) of our method. 
- Missing_quartets_on_trees: This function performs Step (5) of our method. 
- Tallying_quartets_on_trees: This function performs Step (6) of our method. 

