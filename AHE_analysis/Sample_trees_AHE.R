
#required Packages
library(MSCquartets)
library(phylotools)
library('Quartet')
library(tictoc)

#Read files
 
trees=read.tree("desmo900_pisgah_233trees_122individuals.tre")  
Tx=read.table("desmo900_pisgah_233trees_122individuals.csv",sep = ",") 

#Setup to relabel taxa by species
taxa=c()
qE=c()
qD=c()
qC=c()
qG=c()
mC=c()
mE=c()
mG=c()
Outgroup=c()
newtips=function(x)
{
  nt=c()
  for(i in 1:length(x))
  {
    nt[i]=Tx[which(x[i]==Tx[,1]),3]
  }
  newtips=nt
}
for(i in 1:(dim(Tx)[1]))
{
  if(Tx[i,3]=="marm_C"){mC=c(mC,as.character(Tx[i,1]))}
  if(Tx[i,3]=="quad_C"){qC=c(qC,as.character(Tx[i,1]))}
  if(Tx[i,3]=="marm_E"){mE=c(mE,as.character(Tx[i,1]))}
  if(Tx[i,3]=="marm_G"){mG=c(mG,as.character(Tx[i,1]))}
  if(Tx[i,3]=="quad_G"){qG=c(qG,as.character(Tx[i,1]))}
  if(Tx[i,3]=="Outgroup"){Outgroup=c(Outgroup,as.character(Tx[i,1]))}
  if(Tx[i,3]=="quad_D"){qD=c(qD,as.character(Tx[i,1]))}
  if(Tx[i,3]=="quad_E"){qE=c(qE,as.character(Tx[i,1]))}
}
for(i in 1:length(trees))
{
  taxa=c(trees[[i]]$tip.label)
}
taxa=unique(taxa)
qE=intersect(qE,taxa)
qD=intersect(qD,taxa)
qC=intersect(qC,taxa)
qG=intersect(qG,taxa)
mC=intersect(mC,taxa)
mE=intersect(mE,taxa)
mG=intersect(mG,taxa)
Outgroup=intersect(Outgroup,taxa)

alltx=list(qE,mE,mC,qG,qD,qC,mG)  
 

pts=list()
TREES=list()
usedtrees=c()
tis=1
combinations=expand.grid(qE,mE,mC,qG,qD,qC,mG);
dim(combinations)
taxa=c("quad_E", "quad_D", "marm_C", "quad_G", "marm_E","quad_C","marm_G") 



#Sub-sampling procedure
n=10000 #Possible label combinations
r=sample(1:(dim(combinations)[1]),n)#sample
for(i in 1:n)
{ 
  if( (i %% 100)==0 ){cat(i," ")}
  for(j in 1:length(trees))
  {
    tips=intersect(unlist(lapply(combinations[r[i],], as.character)), trees[[j]]$tip.label)
    if(length(tips)==7)
    {
      t=keep.tip(trees[[j]],tips)
      new_tiplabels=newtips(t$tip.label)
      t$tip.label <- new_tiplabels
      TREES[[tis]]=t
      tis=tis+1
    }
  }
}
 
LT=unlist(lapply(TREES, write.tree ))
write.table(LT,paste0("AHE_data.phylo"),quote=F, append = FALSE, sep = " ", dec = ".", row.names = F, col.names = F)
 