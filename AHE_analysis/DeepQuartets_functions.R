### Functions for the method described in the work"
### 'Resolving high-level phylogenetic networks with repeated hybridization in a complex of polytypic salamanders (Caudata: Desmognathus)'

###This R file contain four functions Quartet_Analysis() Trees_quartet_display(), Missing_quartets_on_trees() and Tallying_quartets_on_trees(). These functions correspond to steps 3,4,5,6 from our method respectively. Setps 1 and 2 can be obtained via the function nanuq() of the MSCquartets R package


### Required documentation for the functions
require(MSCquartets)
load("7taxaTrees.RData")
setwd('YOURPATH')
### Function: Quartet_Analysis (Step 3 of hour method)
### Input: (1) a table with the form of the output of the function NANUQ in the package 'MSCquartets' (2) a p-value. The table contains the information about quartet tree concordance factor and the p-value is to termine whether a quarte is deep or non-deep according to the hypothesis test done by the function NANUQ.
Quartet_Analysis=function(pTable,pvalue=.001)
{
  S=pTable
  Qt=list()
  n=dim(S)[1]
  m=dim(S)[2]-1
  NonDeep=c("Trees")
  Deep_hybrid=c("T1","T2")
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
      NonDeep=rbind( NonDeep,Qt[[i]])
    }else
    {
      Deep_hybrid=rbind(Deep_hybrid,Qt[[i]])
    }
    
  }
 
  Quartet_Analysis=list(NonDeep,Deep_hybrid)
}

 
 
 ### Function: Trees_quartet_display (Step 4 of hour method)
 ### Input: The output of the function  Quartet_Analysis
Trees_quartet_display=function(Deep,taxonNames)
{
  Qtree=as.vector(unlist(Deep[[1]]))
  Qnet= unlist(Deep[[2]]) 
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
    
    for(k in 1:length(Qnets[[1]])) 
    {
      for(j in 1:35)
      {
        if(all.equal(unroot(quart[[i]][[j]]),Qnets[[1]][[k]])==T)
        {
          counting=counting+1
        }else if(all.equal( unroot(quart[[i]][[j]]),Qnets[[2]][[k]])==T){
          counting=counting+1
        }
      }
    }
    
    
    cts=c(cts,counting)
    if(counting==35)
    {
      winners=c(i,winners)
    }
    counting=0
  }
  
  T0=c()
  for(i in 1:length(winners))
  {
    T0[i]=write.tree(s7trees[[i]])
  }
  Trees_quartet_display=list(T0,winners,quart,Qnets)
  
}



### Function: Missing_quartets_on_trees (Step 5 of hour method)
### Input: The output of the function  Trees_quartet_display
Missing_quartets_on_trees=function(T0)
{
  QT=list()
  Qnets=T0[[4]]
  quart=T0[[3]]
  possiblewinners=T0[[2]]
 
  powerful=c()
  powerfultree=c()
  cco=0
  otter=list()
  ss=1
  ts=1
  worstcase=c()
  Racks=list()
  winners_Net_list=list()
  for(winwin in  possiblewinners)
  {
    
    Dr=c()
    plow=winwin 
    WQ=quart[[plow]]  
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
    
    for(i in 1:length(Qnet_noTree))
    {
      if(!is.null( Qnet_noTree[[i]])){
      Dr =c(Dr,write.tree( Qnet_noTree[[i]]) )}
      
    }
    Racks[[ts]]= Qnet_noTree
    QT[[ts]]=Dr
    winners_Net_list[[ts]]=winners_Net
    ts=ts+1
    
  }
  
  Missing_quartets_on_trees=list(QT, Racks,winners_Net_list)
  
}
  
### Function:Tallying_quartets_on_trees (Step 6 of hour method) 
### Input: The the first element of the list obtained by the output of the function Missing_quartets_on_trees
Tallying_quartets_on_trees=function(T0,QT,taxonNames)
{
  Qnets=T0[[4]]
  quart=T0[[3]]
  possiblewinners=T0[[2]]
  
  cco=0
  otter=list()
  ss=1
  es=1
  mole=c()
  powerful=c()
  powerfultree=c()
  cco=0
  otter=list()
  ss=1
  ts=1
  winners_Net_list=QT[[3]]
  Racks=QT[[2]]
  worstcase=c()
  for(winwin in  possiblewinners)
  {  
    Qnet_noTree =Racks[[ts]] 
    winners_Net=winners_Net_list[[ts]]
    cco=cco+1
    worstcase[cco]=sum(!is.na(winners_Net))
 
    cts=c()
      for(i in possiblewinners)
      {
        counting=0
        
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
      }
      mole=c(mole,cts)
      otter[[ss]]=cts
      ss=ss+1
      powerfultree=c(powerfultree,winwin)
      powerful=c(powerful,max(cts))

    ts=ts+1
  }
  load("7taxaTrees.RData")
  for(i in 1:945)
  {
    tiplab=s7trees[[i]]$tip.label
    s7trees[[i]]$tip.label=taxonNames[match(c("a","b","c","d","e","f","g"), tiplab)]
  }
  opttrees=c()
  if(length(powerfultree)>0)
  {
    disptrees=powerfultree[powerful==max(powerful)]
    Disptrees=s7trees[disptrees]
    
    for(i in 1:length(disptrees))
    {
      opttrees =c( opttrees,write.tree( s7trees[disptrees[i]] ) ) 
    }
  }else{
    Disptrees="Did not converged"
    opttrees="Did not converged"
  }
  alternativetrees=s7trees[possiblewinners]
  Displayed_Trees=list(Disptrees,alternativetrees)
  Tallying_quartets_on_treess=list( opttrees,Disptrees,alternativetrees)
}
  
 