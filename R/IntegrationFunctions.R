#Integration functions for merging with expression data

mvp.correlate<-function(mvp.tab,exp.tab,method="spearman",fdr.method="BH",verbose=F){
  
  
  #The idea is to have two columns - one the mvp id, the other the gene it is mapped to
  
  
  id.int<-intersect(mvp.tab$gene,exp.tab$gene)
  mvp.tab<-mvp.tab[match(id.int,mvp.tab$gene),]
  exp.tab<-exp.tab[match(id.int,rownames(exp.tab)),]
  
  #Merge
  
  ID.ref<-mvp.tab[,c(1,2)]
  ID.ref<-merge(ID.ref,exp.tab,by.x="gene",by.y=0,all.x=T,sort=F)
  
  #Set up matrices for the correlation function then
  
  mvp.tab<-mvp.tab[,3:ncol(mvp.tab)]
  exp.tab<-ID.ref[,3:ncol(ID.ref)]
  idx.tab<-mvp.tab[,c(1:2)]
  
  #Run correlation function
  
  cx<-correlate.extended(mvp.tab,exp.tab,method=method,verbose=verbose,method.fdr=method.fdr)
  cx<-cbind(idx.tab,cx)
  return(cx)
  
}



#Correlation workhorse function

correlate.extended<-function(obj1,obj2,method,fdr.method,verbose=F){
  
  obj1<-data.matrix(obj1)
  obj2<-data.matrix(obj2)
  estimate<-numeric(); p.value<-numeric()
  
  
  for(i in 1:nrow(obj1)){
    if(verbose){print(i)}
    o1.vec<-as.numeric(as.character(obj1[i,]))
    o2.vec<-as.numeric(as.character(obj2[i,]))
    correlation<-cor.test(o1.vec,o2.vec,method=method)
    p.value[[i]]<-correlation$p.value
    estimate[[i]]<-correlation$estimate
  }
  
  df<-data.frame(estimate=estimate,p.value=p.value,q.value=p.adjust(p.value,method=fdr.method))
  
  return(df)
  
}

#DMR integration function

dmr.correlate<-function(dmr.tab,exp.tab,method,fdr.method,verbose){
  
  dmr.tab2<-dmr.tab
  rownames(dmr.tab)<-strsplit(rownames(dmr.tab),fixed=T,split=".")

  
  
  
  
  
  
  
}


#fDNR collapse function

fdmr.collapse<-function(DMRtab,values,measure,p){
  
  PSum<-DMRtab%>%
    dplyr::select(code,ID)%>%
    mutate(code=as.character(code),
           ID=as.character(ID))
  
  values<-data.frame(values)
  values<-values[rownames(values) %in% PSum$ID,]
  Vals<-merge(PSum,values,by.x="ID",by.y=0,all.x=T)
  dmr.id<-unique(Vals$code)
  Vals<-Vals[,2:ncol(Vals)]
  
    
  if(measure=="mean") {Vals<-Vals%>%
                         group_by(code)%>%
                         summarise_each(funs(mean))
  }
  
  
  if(measure=="median")  {Vals<-Vals%>%
                            group_by(code)%>%
                            summarise_each(funs(median))
  }
  
  
  if(measure=="quantile")  {Vals<-Vals%>%
                              group_by(code)%>%
                              summarise_each(funs(quantile(.,p)))
  }
  
  
  
  Vals<-Vals[,2:ncol(Vals)]%>%
    data.matrix(.)
  
  rownames(Vals)<-dmr.id
  
  return(Vals)
  
  
}
