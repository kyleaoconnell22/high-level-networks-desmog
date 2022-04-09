### Functions for the method described in the work"
### 'Resolving high-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Caudata: Desmognathus)'

###This R file contain two functions Quartet_Analysis() and Displayed_Trees(). This are implementation of the method introduced in the work 'Resolving high-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Caudata: Desmognathus)'
###At the end there is a small example of how these function for a simulated data set


### Required documentation for the functions
require(MSCquartets)
load("7taxaTrees.RData")

### Function: Quartet_Analysis
### Description: This function takes constructs set of quartets identified a deep and non-deep and their associated quartet tree based on quartet gene tree frequencies. Corresponds to Step 3 in the method introduced in the work 'Resolving high-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Caudata: Desmognathus)'
### Input: (1) a table with the form of the output of the function NANUQ in the package 'MSCquartets' (2) a p-value. The table contains the information about quartet tree concordance factor and the p-value is to termine whether a quarte is deep or non-deep according to the hypothesis test done by the function NANUQ.
### Output: A list of two elements: [[1]] the set of quartet trees associated to non-deep quartets [[2]] a matrix containing where each row represent a quartet identified a deep hybrid and 
Quartet_Analysis=function(pTable,pvalue=.001)
{
  S=pTable
  Qt=list()
  n=dim(S)[1]
  m=dim(S)[2]-1
  Qtrees=c("Trees")
  Qnet=c("T1","T2")
  for(i in 1:n)
  { 
    pos=order(as.numeric(S[i,(m-3):(m-1)]),decreasing = T)
    t=pos[1]
    st=pos[2]
    
    if(as.numeric(S[i,m])>pvalue)
    {
      
      if(t==1)
      {
        Qt[[i]]=paste0("((",S[i,1],",",S[i,2],")",",","(",S[i,3],",",S[i,4],"));")
      }else if(t==2)
      {
        Qt[[i]]=paste0("((",S[i,1],",",S[i,3],")",",","(",S[i,2],",",S[i,4],"));")
      }else if(t==3)
      {
        Qt[[i]]=paste0("((",S[i,1],",",S[i,4],")",",","(",S[i,2],",",S[i,3],"));")
      }
    }else{
      if(t==1)
      {
        if(st==2)
        {
          Qt[[i]]=c(paste0("((",S[i,1],",",S[i,2],")",",","(",S[i,3],",",S[i,4],"));"),paste0("((",S[i,1],",",S[i,3],")",",","(",S[i,2],",",S[i,4],"));"))
        }else if(st==3)
        {
          Qt[[i]]=c(paste0("((",S[i,1],",",S[i,2],")",",","(",S[i,3],",",S[i,4],"));"),paste0("((",S[i,1],",",S[i,4],")",",","(",S[i,2],",",S[i,3],"));"))
        }
        
      }else if(t==2)
      {
        if(st==1)
        {
          Qt[[i]]=c(paste0("((",S[i,1],",",S[i,3],")",",","(",S[i,2],",",S[i,4],"));"),paste0("((",S[i,1],",",S[i,2],")",",","(",S[i,3],",",S[i,4],"));")) 
        }else if(st==3)
        {
          Qt[[i]]=c(paste0("((",S[i,1],",",S[i,3],")",",","(",S[i,2],",",S[i,4],"));"),paste0("((",S[i,1],",",S[i,4],")",",","(",S[i,2],",",S[i,3],"));"))
        }
        
      }else if(t==3)
      {
        if(st==1)
        {
          Qt[[i]]=c(paste0("((",S[i,1],",",S[i,4],")",",","(",S[i,2],",",S[i,3],"));"),paste0("((",S[i,1],",",S[i,2],")",",","(",S[i,3],",",S[i,4],"));")) 
        }else if(st==2)
        {
          Qt[[i]]=c(paste0("((",S[i,1],",",S[i,4],")",",","(",S[i,2],",",S[i,3],"));"),paste0("((",S[i,1],",",S[i,3],")",",","(",S[i,2],",",S[i,4],"));"))
        }
      }
    }
    
    if(length(Qt[[i]])==1){
      Qtrees=rbind(Qtrees,Qt[[i]])
    }else
    {
      Qnet=rbind(Qnet,Qt[[i]])
    }
    
  }
  Quartet_Analysis=list(Qtrees,Qnet)
}

### Function: Displayed_Trees
### Description: This function constructs a set of trees that are highly likely displayed in a the network. Corresponds to Step 4, 5, and 6 in the method introduced in the work 'Resolving high-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Caudata: Desmognathus)'
### Input: (1) a table with the form of the output of the function NANUQ in the package 'MSCquartets'. Current implementation works only for a network on 7 taxa (2) a p-value. (3) A vector containing the name of the 7 taxa in the network. This must agree with the pTable.  
### Output: A list of two elements: [[1]] a phylo object containing the tree satisfying all steps of the method, label "Did not converged" if no tree satisfies all steps. [[2]] Teh list of trees satisfying obtained from step 4 of the method. This trees are likely displayed but we are less confident of that.   


Displayed_Trees=function(pTable,pvalue=.05,taxonNames)
{
  Sa=quartetTablePrint(pTable)
  sa=Quartet_Analysis(Sa,pvalue)
  Qtree=as.vector(unlist(sa[[1]]))
  Qnet= unlist(sa[[2]]) 
  Qtree=Qtree[-1]
  Qtrees=list()

  for(i in 1:length(Qtree))
  {
    Qtrees[[i]]=unroot(read.tree(text=Qtree[i]))
  }
  Qnet=Qnet[-1,]
  Qnets=list(list(),list())
  for(i in 1:dim(Qnet)[1])
  {
    Qnets[[1]][[i]]=unroot(read.tree(text=Qnet[i,1]))
    Qnets[[2]][[i]]=unroot(read.tree(text=Qnet[i,2]))
  }
  load("7taxaTrees.RData")
  for(i in 1:945)
  {
    tiplab=s7trees[[i]]$tip.label
    s7trees[[i]]$tip.label=taxonNames[match(c("a","b","c","d","e","f","g"), tiplab)]
  }
  
  quart=list()
  quartets=combn(taxonNames, 3, simplify = FALSE)
  for(i in 1:945)
  {
    t=s7trees[[i]]
    quart[[i]]=list()
    for(j in 1:35)
    {
      quart[[i]][[j]] =drop.tip(t, quartets[[j]] ,rooted=F)
    }
  }

  counting=0
  winners=c()
  cts=c()
  for(i in 1:945)
  {
    for(k in 1:length(Qtrees)) 
    {
      for(j in 1:35)
      {
        if(all.equal( unroot(quart[[i]][[j]]),Qtrees[[k]])==T)
        {
          counting=counting+1
        }
      }
    }
    cts=c(cts,counting)
    if(counting==length(Qtrees))
    {
      winners=c(i,winners)
    }
    counting=0
  }
  
  possiblewinners=winners
  powerful=c()
  powerfultree=c()
  cco=0
  otter=list()
  ss=1
  worstcase=c()
  for(winwin in  possiblewinners)
  {
    plow=winwin 
    WQ=quart[[plow]]  
    # cat(cco/length(possiblewinners),"\n")
    cco=cco+1
    counting_Net=0
    winners_Net=c()
    Qnet_noTree=list() 
    
    for(k in 1:length(Qnets[[1]])) 
    {
      for(j in 1:35)
      {
        if(all.equal( unroot(WQ[[j]]),Qnets[[1]][[k]])==T)
        {
          winners_Net[k]=1
          Qnet_noTree[[k]]=Qnets[[2]][[k]]
        }else if(all.equal( unroot(WQ[[j]]),Qnets[[2]][[k]])==T)
        {
          winners_Net[k]=2
          Qnet_noTree[[k]]=Qnets[[1]][[k]]
        }
      }
    }
    worstcase[cco]=sum(!is.na(winners_Net))
    if(!anyNA(winners_Net))
    { 
      new_things=c()
      for(i in 1:length(winners_Net))
      {
        new_things[i]=Qnet[i,winners_Net[i]]
      }
      counting=0
      winners=c()
      cts=c()
      for(i in 1:945)
      {
        for(k in 1:length(Qnet_noTree)) 
        {
          for(j in 1:35)
          {
            if(all.equal( unroot(quart[[i]][[j]]),Qnet_noTree[[k]])==T)
            {
              counting=counting+1
            }
          }
        }
        cts=c(cts,counting)
        if(counting==length(Qnet_noTree))
        {
          winners=c(i,winners)
        }
        counting=0
      }
      winners
      otter[[ss]]=cts
      ss=ss+1
      which.max(cts)
      sort(cts)
      order(cts)
      powerfultree=c(powerfultree,winwin)
      powerful=c(powerful,max(cts))}
  }
  if(length(powerfultree)>0)
  {
    disptrees=powerfultree[powerful==max(powerful)]
    Disptrees=s7trees[disptrees]
  }else
  {
    Disptrees="Did not converged"
  }
  alternativetrees=s7trees[possiblewinners]
  Displayed_Trees=list(Disptrees,alternativetrees)
}


##### Example

### Simulated data set. Simulated using Hybrid Lambda (Zhu et al 2013)
file="Example_network"  
banyan="(((1:2,h_2#.4:1.5)R:3,((((2:.5)h_2#.4:1,4:1.5)ff:.5,(3:1,(5:.5)h_1#.7:.5)g:1)e:1,((h_1#.7:.5, 6:1)h:1,7:2)f:1)c:2)a:2)r;"
plot(read.evonet(text=banyan),arrows = F,use.edge.length = TRUE)

## Run NANUQ (Steps 1 and 2)
genetrees <- read.tree(file)
pTable=NANUQ(genetrees, outfile = file,alpha=.0001, beta=.95)
 
## Run Displayed_Trees Steps 3, 4, 5, and 6
Trees=Displayed_Trees(pTable,pvalue=.0001,taxonNames=c("1_1","2_1","3_1","4_1","5_1","6_1","7_1"))

 
