### Required documentation for the procedure
require(MSCquartets)
library(phylotools)
library('Quartet')

load("7taxaTrees.RData") #Rfile containig all tree topologies on 7 taxa
#Open this file and setwd in the functions file
source("DeepQuartets_functions.R", echo=F)
#We show the performance of the procedure by looking at simulated gene trees from the following network
network="(((1:2,h_2#.4:1.5)R:3,((((2:.5)h_2#.4:1,4:1.5)ff:.5,(3:1,(5:.5)h_1#.7:.5)g:1)e:1,((h_1#.7:.5, 6:1)h:1,7:2)f:1)c:2)a:2)r;"
plot(read.evonet(text=network),arrows = F,use.edge.length = TRUE)

# Simulated data set. Simulated using Hybrid Lambda (Zhu et al 2013) from the Network above
file="Example_network_sample"  
genetrees <- read.tree(file)

# We now show each of the steps of the procedure to obtain the displayed trees:

## (Steps 1 and 2) 
####Run NANUQ on simulated set of genes 
pTable=NANUQ(genetrees, outfile = file,alpha=.0001, beta=.95)
Sa=quartetTablePrint(pTable)
## (Step 3) 
####The following function takes the output from Step 2 and classifies all quartets as wither non-deep or deep-hybrid and shows the tree (or trees) associated to such quartet
Deep=Quartet_Analysis(Sa,pvalue=.0001);
Deep
## (Step 4)
#### The following function will construct all trees that display all trees associated with non-deep quartet networks, and one of the quartet trees associated to each quartet identified as deep-hybrid.
T0=Trees_quartet_display(Deep,taxonNames=c("1_1","2_1","3_1","4_1","5_1","6_1","7_1"))
T0[[1]]
## (Step 5)
#### 
QT=Missing_quartets_on_trees(T0)
QT[[1]]
## (Step 6)
#### 
ST=Tallying_quartets_on_trees(T0,QT,taxonNames=c("1_1","2_1","3_1","4_1","5_1","6_1","7_1"))
ST[[1]]
for(i in 1:length(ST[[1]]))
{
  plot(ST[[2]][i])
}



