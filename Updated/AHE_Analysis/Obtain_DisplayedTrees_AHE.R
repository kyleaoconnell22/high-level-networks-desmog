### This file contains the procedure to obtain the Backbone tree for the network in the work
### 'Resolving high-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Caudata: Desmognathus)'

### Required documentation for the procedure
require(MSCquartets)
load("7taxaTrees.RData") #Rfile containig all tree topologies on 7 taxa
source("DeepQuartets_functions.R", echo=F) 
library(MSCquartets)
library(phylotools)
library('Quartet')
library(tictoc)


#load trees obtained in file 'Create_subsmaple_for_nanu_AHE.R' 
trees=read.tree("AHE_data.phylo")#AHE 381
 
 
## (Steps 1) 
####We use the Function 'quartetTableParallel()' which is the procedure in NANUQ that does the quartet tallying (step 1) 
taxa=c("quad_E", "quad_D", "marm_C", "quad_G", "marm_E","quad_C","marm_G") 
pTable  = quartetTableParallel( trees,taxonnames =taxa,numCores=15)
pTable  = quartetTableResolved(pTable) # we removed unresolved as there were no unresolved quartets

## (Steps 2) 
####We use the Function 'quartetTreeTestInd()' which is the procedure in NANUQ that applies a statistical hypothesis test to the empirical concordance factors  (step 2)
pTable  = quartetTreeTestInd(pTable , "T3")
pTable  = quartetStarTestInd (pTable ) # We test for unresolved quartets to be cautios
save(pTable , file=paste0("pTable_AHE.RData"))
Sa=quartetTablePrint(pTable)

## (Step 3) 
####The following function takes the output from Step 2 and classifies all quartets as wither non-deep or deep-hybrid and shows the tree (or trees) associated to such quartet
Deep=Quartet_Analysis(Sa, pvalue=.005);
Deep
save(Deep , file=paste0("Deep_AHE.RData"))
## (Step 4)
#### The following function will construct all trees that display all trees associated with non-deep quartet networks, and one of the quartet trees associated to each quartet identified as deep-hybrid.
T0=Trees_quartet_display(Deep,taxonNames= c("quad_E", "quad_D", "marm_C", "quad_G", "marm_E","quad_C","marm_G") )
T0[[1]]
save(T0 , file=paste0("T0_AHE.RData"))
## (Step 5)
#### 
QT=Missing_quartets_on_trees(T0)
QT[[1]]
 
## (Step 6)
#### 
DisplayeTrees=Tallying_quartets_on_trees(T0,QT,taxonNames= c("quad_E", "quad_D", "marm_C", "quad_G", "marm_E","quad_C","marm_G") )
DisplayeTrees[[1]]
for(i in 1:length(DisplayeTrees[[1]]))
{
  plot(ST[[2]][i])
}
 
save(DisplayeTrees, file=paste0("DisplayeTrees_AHE.RData"))



















