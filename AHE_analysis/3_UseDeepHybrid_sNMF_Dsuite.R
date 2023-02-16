###  
### 'Resolving high-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Caudata: Desmognathus)'

taxonNames= c("quad_E", "quad_D", "marm_C", "quad_G", "marm_E","quad_C","marm_G")
#load backbone trees
load("DisplayeTrees_AHE.RData")
BackboneTrees=DisplayeTrees[[2]]
#load information of Deep hybrids
load("T0_AHE.RData")
Deep_hybrid=T0[[4]]


Missing_deep=c()
for(i in 1:length(Deep_hybrid[[1]]))
{
  tx=Deep_hybrid[[1]][[i]]$tip.label
  T4=drop.tip(BackboneTrees[[1]],setdiff(taxonNames,tx))
  
  if(all.equal( unroot(T4),Deep_hybrid[[1]][[i]])==T)
  {
    Missing_deep=c(Missing_deep, write.tree( Deep_hybrid[[1]][[i]]) )
  }else if(all.equal(unroot(T4),Deep_hybrid[[2]][[i]])==T){
    Missing_deep=c(Missing_deep, write.tree( Deep_hybrid[[2]][[i]]) )
  }
  
}

#list of quartets not displayed in the backbone tree:
Missing_deep
 
 
################################
### Adding extra information ###
################################
 
#The admixture signal from sNMF using the AHE data (Pyron et al. 2022c) suggests that there is a hybrid event between Desmognathus marmoratus and "quad" G. 

### This agrees with quartets  [19],and [20]
Missing_deep[19]
Missing_deep[20]
#The sNMF analysis of the  GBS data (Fig. 2) shows that "quad" G is a mixture of D. marmoratus, D. mavrokoilius, and the group that represents D. intermedius + "quad" C. This is the scenario depicted by edges a and b

 ### This agrees with quartets [16,26,27].
Missing_deep[16]
Missing_deep[26]
Missing_deep[27]
 
#From the same admixture analyses, we obtain two other hybrid events. One suggests that "marm" G is a mixture of Desmognathus intermedius and D. marmoratus
 
### This agrees with quartet  [12] 
Missing_deep[12]

## when combined with the other events already mentioned, agrees with quartets [13,23,25]
Missing_deep[13]
Missing_deep[23]
Missing_deep[25]

#The final hybridization event obtained from sNMF analysis of the GBS data suggests that "quad" C is a mixture of Desmognathus intermedius and small amounts of D. marmoratus and D. mavrokoilius

### This agrees with quartets  [1,9,10,17,18]
Missing_deep[1]
Missing_deep[9]
Missing_deep[10] 
Missing_deep[17]
Missing_deep[18]

# We suspect the signal of D. mavrokoilius could come from the hybridization between D. mavrokoilius and D. marmoratus

#When combining this hybridization and the other described above such a network has 4 more rows [11,14,15,24]
Missing_deep[11]
Missing_deep[14]
Missing_deep[15] 
Missing_deep[24]


