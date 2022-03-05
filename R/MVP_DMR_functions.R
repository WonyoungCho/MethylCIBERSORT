Annotate.object <- function(object, use.dB = "mean"){

  data("probe_450K_VCs_af.RData")
  data("probe_450k_features_V4.rda")
  data("Epic_450k_consensus.RData")
  data("450k_whitelist.RData")
  probeData <- probe.features


  #compute dB

  if(use.dB == "mean") { object$dB <- object[,9] - object[,8] ; message(paste0("Computing diff"," ", colnames(object)[[9]] ," - " ,colnames(object)[[8]]))}
  if(use.dB == "median") { object$dB <- object[,11] - object[,10] ; message(paste0("Computing diff"," ", colnames(object)[[11]] ," - " ,colnames(object)[[10]]))}

  #data(probeData)

  if(class(object)=="list") {

    object <- lapply(object, function(x) merge(x,probeData,by = "ID"))

  } else {

    object <- merge(object, probeData, by="ID")

  }

  return(object)

}


find.mvp<-function(values, type=c("beta","M"), design, coef, contrast.matrix,classes, measure.type,TSBH=T,alpha.TSBH=0.05,... ){


  #Create results list

  Results.List <- list()

  #This block of code generates an object called Fit using limma

  message("Fitting Model")
  if(is.null(contrast.matrix)){

    Fit<-lmFit(values,design)%>%
      eBayes(.,trend=T)


  } else {

    Fit<-lmFit(values,design)%>%
      contrasts.fit(.,contrast.matrix)%>%
      eBayes(.,trend=T)

  }

  # Add to list here

  Results.List$Fit <- Fit



  if(length(coef) == 1) {

    tt <- topTable(Fit, coef= coef, number = nrow(values))%>%transform(.,ID=rownames(.))

    if(type=="beta"){

      if(is.null(classes)){stop("computing delta beta requires a vector of groups")}
      message("Computing beta value summary statistics ")

      dB.mean <-sumClass(classes=classes,values=values,measure="mean")%>%
        data.frame(.)

      colnames(dB.mean) <- paste0("mean",".",colnames(dB.mean))

      dB.median <- sumClass(class=classes, values=values, measure="median")%>%
        data.frame(.)

      colnames(dB.median) <- paste0("median",".", colnames(dB.median))

      dB.mean$ID <- rownames(dB.mean) ; dB.median$ID <- rownames(dB.median)

      dB.tab <- merge(dB.mean,dB.median,by="ID")


      tt<-merge(tt,dB.tab,by="ID")

    }



    if(TSBH) {  tt <- TSBH.adjust(topTable=tt, alpha = alpha.TSBH) }
    Results.List$tt <- tt

  } else {


    ## This block of code is designed to work with multiple specified coefficients
    message("Multiple coefficients specified, preparing limma topTables")

    topTables.list <- list()


    for( i in 1:length(coef)) {

      tt <- topTable(Fit, coef= i,number = nrow(values))%>%transform(., ID=rownames(.))

      if(type=="beta") {

        if(is.null(classes)){stop("computing delta beta requires a vector of groups")}
        message("Computing beta value summary statistics ")

        cl.v <- classes[[i]]
        dB.mean <-sumClass(classes=cl.v,values=values,measure="mean")%>%
          data.frame(.)

        colnames(dB.mean) <- paste0("mean",colnames(dB.mean))

        dB.median <- sumClass(classes = classes[[i]], values=values, measure="median")%>%
          data.frame(.)

        colnames(dB.median) <- paste0("median", colnames(dB.median))

        dB.tab <- merge(dB.mean,dB.median,by=0)%>%transform(.,ID=rownames(.))

        tt <-merge(tt,dB.tab,by="ID")

      }

      topTables.list[[i]] <- tt
      message(i)

    }

    topTables.list <- lapply(topTables.list, function(x) transform(x, ID=rownames(x)))

    if(TSBH) {
      for(i in 1:length(topTables.list)) { topTables.list[[i]] <- TSBH.adjust(topTable=tt[[i]],alpha=alpha.TSBH)}
    }



    Results.List$tt <- topTables.list


  }



  return(Results.List)
}



# New VMP stuff



find.vmp<-function(values, type=c("beta","M"), design, coef, contrast.matrix,TSBH=T,alpha.TSBH=0.05,... ){

  if(type=="beta"){ values = 2^(log2(values)-log2(1-values)); message("Converting to M values for variance modelling")}

  Results.List <- list()

  #This block of code generates an object called Fit using limma

  if(is.null(contrast.matrix)){

    Fit<-varFit(values,design)%>%
      eBayes(.)



  } else {

    Fit<-varFit(values,design)%>%
      contrasts.fit(.,contrast.matrix)%>%
      eBayes(.,trend=T)


  }

  Results.List$Fit <- Fit


  if(length(coef) == 1 ) {

    tt <- topVar(Fit,number=nrow(values),coef=coef)%>%
      transform(.,ID=rownames(.))

    if(TSBH) { tt <- TSBH.adjust.v(tt, alpha = alpha.TSBH)}

    Results.List$tt <- tt

  } else {

    tt.list <- list()

    for ( i in 1:length(coef)) {

      tt.list[[i]] <- topVar(Fit, number= nrow(values), coef=coef[[i]])%>%
        transform(., ID=rownames(.))

    }


    if(TSBH) {

      for( i in 1:length(coef)) { tt.list[[i]] <- TSBH.adjust.v(topTable= tt.list[[i]],alpha=alpha.TSBH)}

    }

    Results.List$tt <- tt.list
  }


  return(Results.List)

}





# This is the class summarisation function
sumClass<-function(classes,measure=c("median","mean"),values,verbose=T){

  levs<-levels(factor(classes))
  output<-list()

  for(i in 1:length(levs)){

    class<-classes %in% levs[[i]]
    val<-values[,class]

    if(measure=="median"){output[[i]]<-rowMedians(val)}
    if(measure=="mean")  {output[[i]]<-rowMeans(val)}
  }

  df<-do.call(cbind,output)
  colnames(df)<-levs
  rownames(df)<-rownames(values)
  return(df)

}

# This is the TSBH adjustment function
TSBH.adjust<-function(topTable,alpha){

  Fit<-topTable
  qvals<-mt.rawp2adjp(rawp=Fit$P.Value,proc="TSBH",alpha)
  q<-qvals$adjp[,2]
  index<-qvals$index
  Fit$adj.P.Val<-q[order(index)]
  return(Fit)
}

TSBH.adjust.v <- function(topTable,alpha){

  Fit<-topTable
  qvals<-mt.rawp2adjp(rawp=Fit$P.Value,proc="TSBH",alpha)
  q<-qvals$adjp[,2]
  index<-qvals$index
  Fit$Adj.P.Value<-q[order(index)]
  return(Fit)
}


#This implements the fDMR function

fDMR<-function(signif.only=F,mvp.tab,mvp.fdr,n.mvp,dmr.fdr,mvp.dB,dmr.dB,beta,dB.measure,classes,exclude.regions,os.p=F,resolution=c("high","low")){

  #Initial analysis...

  #Here, I aim to reduce stuff to candidate sigDMRs based on cutoff criteria

  #Annotate first with feature group
  message("annotating with feature groups")

  if(resolution=="low") {

  mvp.tab<-mvp.tab%>%
    mutate(feature2=as.character(feature),feature2=ifelse(feature2 %in% c("TSS1500","TSS200","1stExon"),
                                                          "promoter",ifelse(feature2 %in% "Body",
                                                                            "Body",feature2)),
           code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(adj.P.Val))

  } else {


    mvp.tab <- mvp.tab%>%mutate(feature2=as.character(feature), code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(adj.P.Val))

  }



  if(!is.null(exclude.regions)){

    mvp.tab<-filter(mvp.tab,!feature2 %in% exclude.regions)

  }


  mvp.cutoff<-mvp.tab%>%filter(adj.P.Val < mvp.fdr, dB > mvp.dB | dB < -mvp.dB)

  #Then carry out fDMR filtration to select dmrs that meet criteria for candidates

  mvp.count<-dplyr::count(mvp.cutoff,code)%>%
    filter(n > (n.mvp-1))


  #Then we reduce the original table to candidate dmrs for dmr calling and stuff
  message("setting up probe table for candidate dmrs")
  if(signif.only==FALSE){mvp.tab2<-mvp.tab%>%filter(code %in% mvp.count$code)}else
  {mvp.tab2<-mvp.cutoff%>%filter(code %in% mvp.count$code)}

  mvp.tab2<-mvp.tab2[order(mvp.tab2$code),]
  beta.tab2<-beta[match(mvp.tab2$ID,rownames(beta)),]


  if(os.p==T){

    #One sided P value conversion here...

    message("Converting p values to 1 sided values for p.value combination")
    P.os<-mapply(function(t,p) p.ts.os(t,p),mvp.tab2$t,as.character(mvp.tab2$P.Value))
    mvp.tab2$P.os<-P.os
  }  else {

    message("parameter choice is to not conver to one-sided p values")
    mvp.tab2$P.os<-as.numeric(as.character(mvp.tab2$P.Value))

  }

  #Then do the whole p value combination shenanigans.

  message("SLK correction")
  dmr.beta <- split(data.frame(beta.tab2), factor(mvp.tab2$code))
  corel <- lapply(dmr.beta, function(x) cor(t(x)))
  weights <- lapply(corel, function(x) 1/apply(x^2,1,sum))

  if(os.p==T){
    dmr.ind.p <- split(mpfr(mvp.tab2$P.os), factor(mvp.tab2$code))
    dmr.ind.p<-lapply(dmr.ind.p,function(x) asNumeric(x))
    dmr.qp <- lapply(dmr.ind.p, qnorm)
    dmr.qp.w <- mapply("*", dmr.qp, weights)
  } else {

    dmr.ind.p<- split(mvp.tab2$P.os, factor(mvp.tab2$code))
    dmr.qp <- lapply(dmr.ind.p, qnorm)
    dmr.qp.w <- mapply("*", dmr.qp, weights)

  }

  #Then generate p value distributions

  if(class(dmr.qp.w) == "matrix")
  {
    dmr.stat <- sum(dmr.qp.w)
  }else
  {
    dmr.stat <- lapply(dmr.qp.w, sum)
  }

  message("calculating DMR p values")
  dmr.sd <- lapply(weights, function(x) sqrt(sum(x^2)))
  dmr.p <- mapply(function(x,y) pnorm(x,0, sd=y), dmr.stat, dmr.sd)
  if(os.p==T) {dmr.p <- lapply(as.character(dmr.p), function(x) p.os.ts(x))}
  dmr.p <- unlist(dmr.p)

  #Once this is done, I can go on and build up a table

  code.p<-unique(factor(mvp.tab2$code))
  dmr.summary<-data.frame(code=code.p,dmr.p=dmr.p,dmr.FDR=p.adjust(dmr.p,method="BH"))

  message("Exporting results tables")
  mvp.summary<-merge(dmr.summary,mvp.tab2,by="code")

  #Then further annotate results table.

  dmr.summary<- filter(dmr.summary,dmr.FDR < dmr.fdr)
  mvp.summary<- filter(mvp.summary,code %in% dmr.summary$code, adj.P.Val < mvp.fdr, dB > mvp.dB | dB < -mvp.dB )

  dB.summary<-mvp.summary%>%
    group_by(code)%>%
    summarise(Mean.dB=mean(dB),Median.dB=median(dB))

  dmr.summary<-merge(dmr.summary,dB.summary,by="code",all.x=T)
  if(dB.measure=="mean"){dmr.summary<-filter(dmr.summary,Mean.dB > dmr.dB | Mean.dB < -dmr.dB)}
  if(dB.measure=="median"){dmr.summary<-filter(dmr.summary,Median.dB > dmr.dB | Median.dB < -dmr.dB)}


  mvp.summary<- filter(mvp.summary,code %in% dmr.summary$code, adj.P.Val < mvp.fdr, dB > mvp.dB | dB < -mvp.dB )


  return(list(Probes=mvp.summary,dmrs=dmr.summary))

}

#Coxfit.MVP - this function basically estimates cox model statistics for matrices of beta or m values.

coxfit.mvp<-function(time,status,strata.f,values,weighted=FALSE,verbose=T){

  values<-t(values)
  sData<-Surv(time,status)


  if(weighted==FALSE){

    results.list<-list()


    for(i in 1:ncol(values)){tryCatch({

      if(is.null(strata.f)){ cph.obj<-coxph(sData~values[,i])
      }  else {

        cph.obj<-coxph(sData~values[,i]+strata(strata.f),singular.ok=T)

      }

      results.list[[i]]<-coef(summary(cph.obj))
      rownames(results.list[[i]])<-colnames(values)[[i]]
      if(verbose==T){print(i)}
    },error=function(cond){})

    }


    l2<-do.call(rbind,results.list)%>%
      data.frame(.)

    colnames(l2)<-c("Coef","HR","se.Coef","Z","p.value")
    l2$ID<-rownames(l2)
    return(l2)

  } else {


    p.values<-c(rep(0,ncol(values)))
    coefs<-c(rep(0,ncol(values)))

    if(is.null(strata)){

      for(i in 1: ncol(values)) {

        fit<-coxphw(sData~values[,i])
        p.values[[i]]<-fit$prob
        coefs[[i]]<-fit$coef
      }
    } else {

      for(i in 1: ncol(values)){
        fit<-coxphw(sData~values[,i]+strata.f)
        p.values[[i]]<-fit$prob
        coefs[[i]]<-fit$coef
      }
    }

    ID.vec<-colnames(values)
    df<-data.frame(ID=ID.vec,coef=coefs,p.value=p.values)%>%
      mutate(HR=exp(coef))
    return(df)

  }
}

#Develop a version of fDMR for survival

surv.fDMR<-function(beta,mvp.cox,alpha=0.05,TSBH=T,mvp.fdr=0.05,dmr.fdr,p.os,signif.only=F,n.mvp=3,exclude.regions,resolution=c("high","low")){

  #This first bit performs TSBH correction if required

  if(TSBH){

    Fit<-mvp.cox
    qvals<-mt.rawp2adjp(rawp=Fit$p.value,proc="TSBH",alpha)
    q<-qvals$adjp[,2]
    index<-qvals$index
    Fit$adj.P.Val<-q[order(index)]

  } else {

    Fit$adj.P.Val<- p.adjust(Fit$adj.P.Val,method="BH")

  }

  mvp.tab<-Fit
  rm(Fit)

  #From the next bit we begin to pass variables into a P.value combination function after filtering.

  #Annotation and Filtration module

  message("annotating with feature groups")
  mvp.tab$ID<-rownames(mvp.tab)
  mvp.tab<-merge(mvp.tab,probe.features,by.x="ID",by.y=0,sort=F,all.x=T)

  if(resolution=="low") {

    mvp.tab<-mvp.tab%>%
      mutate(feature2=as.character(feature),feature2=ifelse(feature2 %in% c("TSS1500","TSS200","1stExon"),
                                                            "promoter",ifelse(feature2 %in% "Body",
                                                                              "Body",feature2)),
             code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(adj.P.Val))

  } else {


    mvp.tab <- mvp.tab%>%mutate(feature2=as.character(feature), code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(adj.P.Val))

  }



  if(!is.null(exclude.regions)){

    mvp.tab<-filter(mvp.tab,!feature2 %in% exclude.regions)

  }



  #The next bit filters by FDR and number of mvp cutoffs

  mvp.cutoff<-mvp.tab%>%filter(adj.P.Val < mvp.fdr)
  mvp.count<-dplyr::count(mvp.cutoff,code)%>%
    filter(n > (n.mvp-1))

  #Then we reduce the original table to candidate dmrs for dmr calling and stuff

  message("setting up probe table for candidate dmrs")
  if(signif.only==FALSE){mvp.tab2<-mvp.tab%>%filter(code %in% mvp.count$code)}else
  {mvp.tab2<-mvp.cutoff%>%filter(code %in% mvp.count$code)}

  mvp.tab2<-mvp.tab2[order(mvp.tab2$code),]
  beta.tab2<-beta[match(mvp.tab2$ID,rownames(beta)),]



  #P value conversion to one sided if needed

  if(os.p==T){

    #One sided P value conversion here...

    message("Converting p values to 1 sided values for p.value combination")
    P.os<-mapply(function(t,p) p.ts.os(t,p),mvp.tab2$t,as.character(mvp.tab2$P.Value))
    mvp.tab2$P.os<-P.os
  }  else {

    message("parameter choice is to not convert to one-sided p values")
    mvp.tab2$P.os<-as.numeric(as.character(mvp.tab2$P.Value))

  }

  #Correct p values for combination purposes
  message("Stouffer Liptak Correction")
  dmr.beta <- split(data.frame(beta.tab2), factor(mvp.tab2$code))
  corel <- lapply(dmr.beta, function(x) cor(t(x)))
  weights <- lapply(corel, function(x) 1/apply(x^2,1,sum))

  if(os.p==T){
    dmr.ind.p <- split(mpfr(mvp.tab2$P.os), factor(mvp.tab2$code))
    dmr.ind.p<-lapply(dmr.ind.p,function(x) asNumeric(x))
    dmr.qp <- lapply(dmr.ind.p, qnorm)
    dmr.qp.w <- mapply("*", dmr.qp, weights)
  } else {

    dmr.ind.p<- split(mvp.tab2$P.os, factor(mvp.tab2$code))
    dmr.qp <- lapply(dmr.ind.p, qnorm)
    dmr.qp.w <- mapply("*", dmr.qp, weights)

  }

  #Then generate p value distributions

  if(class(dmr.qp.w) == "matrix")
  {
    dmr.stat <- sum(dmr.qp.w)
  }else
  {
    dmr.stat <- lapply(dmr.qp.w, sum)
  }

  message("calculating DMR p values")
  dmr.sd <- lapply(weights, function(x) sqrt(sum(x^2)))
  dmr.p <- mapply(function(x,y) pnorm(x,0, sd=y), dmr.stat, dmr.sd)
 if(os.p==T) {dmr.p <- lapply(as.character(dmr.p), function(x) p.os.ts(x))}
  dmr.p <- unlist(dmr.p)

  #Once this is done, I can go on and build up a table

  message("Generating output objects")
  code.p<-unique(factor(mvp.tab2$code))
  dmr.summary<-data.frame(code=code.p,dmr.p=dmr.p,dmr.FDR=p.adjust(dmr.p,method="BH"))
  mvp.summary<-merge(dmr.summary,mvp.tab2,by="code")

  return(list(Probes=mvp.summary,dmrs=dmr.summary))
}


#Function to handle mpfr output

o.mpfr<-function(n) capture.output(n)[2]%>%substr(.,5,nchar(.))

#Functions for 1-sided / 2-sided p value conversions

p.ts.os<-function(t,p){
  p<-mpfr(p,base=10)
  if( t < 0 ) {px= 0.5* (1-p)} else
  { px = p*0.5}

  px<-asNumeric(px)
  return(as.character(px))
}

p.os.ts<-function(p){
  p<-mpfr(p)
  if( p <= 1/2) { p = 2*p} else { p = 2(1-p)}
  p<-asNumeric(p)
  return(p)
}


# fVMR for coordinate variability calling


fVMR<-function(signif.only=F,vmp.tab,vmp.fdr,n.vmp,vmr.fdr,Mtab,exclude.regions,os.p=F,resolution=c("high","low")){

  #Initial analysis...

  #Here, I aim to reduce stuff to candidate sigvmrs based on cutoff criteria

  #Annotate first with feature group
  message("annotating with feature groups")
  if(!"ID" %in% colnames(vmp.tab) ){vmp.tab$ID<-rownames(vmp.tab)}


  if(resolution=="low") {

    vmp.tab<-vmp.tab%>%
      mutate(feature2=as.character(feature),feature2=ifelse(feature2 %in% c("TSS1500","TSS200","1stExon"),
                                                            "promoter",ifelse(feature2 %in% "Body",
                                                                              "Body",feature2)),
             code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(Adj.P.Value))

  } else {


    vmp.tab <- vmp.tab%>%mutate(feature2=as.character(feature), code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(Adj.P.Value))%>%
      dplyr::filter(!feature2 %in% exclude.regions)

  }






  if(!is.null(exclude.regions)){

    vmp.tab<-filter(vmp.tab,!feature2 %in% exclude.regions)

  }


  vmp.cutoff<-vmp.tab%>%filter(Adj.P.Value < vmp.fdr)

  #Then carry out fvmr filtration to select vmrs that meet criteria for candidates

  vmp.count<-dplyr::count(vmp.cutoff,code)%>%
    filter(n > (n.vmp-1))


  #Then we reduce the original table to candidate vmrs for vmr calling and stuff
  message("setting up probe table for candidate vmrs")
  if(signif.only==FALSE){vmp.tab2<-vmp.tab%>%filter(code %in% vmp.count$code)}else
  {vmp.tab2<-vmp.cutoff%>%filter(code %in% vmp.count$code)}

  vmp.tab2<-vmp.tab2[order(vmp.tab2$code),]
  Mtab.tab2<-Mtab[match(vmp.tab2$ID,rownames(Mtab)),]


  if(os.p==T){

    #One sided P value conversion here...

    message("Converting p values to 1 sided values for p.value combination")
    P.os<-mapply(function(t,p) p.ts.os(t,p),vmp.tab2$t,as.character(vmp.tab2$P.Value))
    vmp.tab2$P.os<-P.os
  }  else {

    message("parameter choice is to not convert to one-sided p values")
    vmp.tab2$P.os<-as.numeric(as.character(vmp.tab2$P.Value))

  }

  #Then do the whole p value combination shenanigans.

  message("Stouffer Liptak correction")
  vmr.Mtab <- split(data.frame(Mtab.tab2), factor(vmp.tab2$code))
  corel <- lapply(vmr.Mtab, function(x) cor(t(x)))
  weights <- lapply(corel, function(x) 1/apply(x^2,1,sum))

  if(os.p==T){
    vmr.ind.p <- split(mpfr(vmp.tab2$P.os), factor(vmp.tab2$code))
    vmr.ind.p<-lapply(vmr.ind.p,function(x) asNumeric(x))
    vmr.qp <- lapply(vmr.ind.p, qnorm)
    vmr.qp.w <- mapply("*", vmr.qp, weights)
  } else {

    vmp.tab2$P.os <- as.numeric(vmp.tab2$P.os)
    vmr.ind.p<- split(vmp.tab2$P.os, factor(vmp.tab2$code))
    vmr.qp <- lapply(vmr.ind.p, qnorm)
    vmr.qp.w <- mapply("*", vmr.qp, weights)

  }

  #Then generate p value distributions

  if(class(vmr.qp.w) == "matrix")
  {
    vmr.stat <- sum(vmr.qp.w)
  }else
  {
    vmr.stat <- lapply(vmr.qp.w, sum)
  }

  message("calculating vmr p values")
  vmr.sd <- lapply(weights, function(x) sqrt(sum(x^2)))
  vmr.p <- mapply(function(x,y) pnorm(x,0, sd=y), vmr.stat, vmr.sd)
  if(os.p==T) {vmr.p <- lapply(as.character(vmr.p), function(x) p.os.ts(x))}
  vmr.p <- unlist(vmr.p)

  #Once this is done, I can go on and build up a table

  code.p<-unique(factor(vmp.tab2$code))
  vmr.summary<-data.frame(code=code.p,vmr.p=vmr.p,vmr.FDR=p.adjust(vmr.p,method="BH"))

  message("Exporting results tables")
  vmp.summary<-merge(vmr.summary,vmp.tab2,by="code")

  #Then further annotate results table.

  vmr.summary<- filter(vmr.summary,vmr.FDR < vmr.fdr)
  vmp.summary<- filter(vmp.summary,code %in% vmr.summary$code, Adj.P.Value < vmp.fdr)

  var.summary<-vmp.summary%>%
    group_by(code)%>%
    summarise(Mean.logFC=mean(LogVarRatio),Median.logFC=median(LogVarRatio))

  vmr.summary<-merge(vmr.summary,var.summary,by="code",all.x=T)


  vmp.summary<- filter(vmp.summary,code %in% vmr.summary$code, Adj.P.Value < vmp.fdr)


  return(list(Probes=vmp.summary,vmrs=vmr.summary))

}

