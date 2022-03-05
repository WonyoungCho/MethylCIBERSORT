# Featureselect

# This function simply selects features for a given cell type to permit decon
#Algorithm - combine, normalise, remove blacklist, set up phenotype table, apply limma, apply filter, return.

#' @export

FeatureSelect <- function(RGSet.Cancer, do.BMIQ = F, norm.method = "Funnorm",AllGap = 0.3, nnGap = 0.1, whitelist = NULL, verbose = T, feature.number = 50 , probeFDR = 0.01) {

  #Load libraries and stuff

  require(plyr)
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  require(Biobase)
  require(minfi)
  require(BiocGenerics)
  require(limma)

  #Load RGsets for infiltrating cell line stuff

  data("RGSets_infiltrate.RData")

 if(verbose){ message("Beginning combination and evaluation processes") }

  if(!class(RGSet.Cancer)=="RGChannelSet"){
    stop("Your set of cancer cell lines needs to be an RGChannelSet")
  }

  #Combine and normalise

  Combined.RGSet <- BiocGenerics::combine(RGSet.Cancer,MicSet2)


  if(norm.method == "Funnorm") {

    if(verbose){ message ("Functional normalisation in process")}
    Combined <- preprocessFunnorm(Combined.RGSet, bgCorr = F, dyeCorr = F, verbose = verbose)

  }


  if(norm.method == "SWAN") {

    Combined <- preprocessSWAN(Combined.RGSet)

  }


  if(norm.method == "Raw") {

    Combined <- preprocessRaw(Combined.RGSet)

  }


  #Continue
  Pheno.cancer <- data.frame(Sample = colnames(RGSet.Cancer))%>%
    mutate(Group = "Cancer")

  Pheno <- rbind(Pheno.cancer, Ref.pheno2)
  Beta <- getBeta(Combined)
  rm(Combined.RGSet, RGSet.Cancer, Combined)
  Beta <- Beta[complete.cases(Beta),]

#This module executes BMIQ normalisation

if(!is.null(whitelist)) { Beta <- Beta[rownames(Beta) %in% whitelist ,]}

if(do.BMIQ) {
  message("Performing BMIQ normalisation")
  Beta <- HMkit.BMIQ(vals = Beta, nfit = 10000, plot = F)
}

#From here I can work on actually defining signatures

  design <- model.matrix(~factor(Pheno$Group))
  colnames(design) <- levels(factor(Pheno$Group))
  GroupsVec <- levels(factor(Pheno$Group))

  #This model sweeps through for one vs all comparisons

  if(verbose) { message ("Beginning one versus all sweeps")}

  ModelOutputs.list <- list()

  for(i in 1:length(GroupsVec)) {

    message(paste0("Initiating"," ",i))
    Vec1 <- as.character(Pheno$Group)
    Vec1 <- ifelse(Vec1 %in% GroupsVec[[i]], "One","Other")
    design <- model.matrix(~0+factor(Vec1))
    colnames(design) <- levels(factor(Vec1))

    Fit <- find.mvp(values= Beta, type = "beta", design, coef = 1, contrast.matrix = makeContrasts(One-Other,levels = design),classes = factor(Vec1), TSBH = T, alpha.TSBH = 0.05, measure.type = "median")

    message(paste0("Finishing"," ", i))
    ModelOutputs.list[[i]] <- Fit

  }


  #Filter by initial criteria here

  if(verbose) { message("Filtering probes for further application here")}

  Models2 <- lapply(ModelOutputs.list, function(x) Annotate.object(object = x$tt,use.dB = "median"))
  names(Models2) <- GroupsVec

  Filtered.tt <- lapply(Models2, function(x) filter(x, adj.P.Val < probeFDR, dB > AllGap | dB < -AllGap))

  lapply(Filtered.tt,dim)

  #Then we filter by nnGap - this ensures there is some degree of separation

  Medians.List <- list()

  for (i in 1:length(GroupsVec)) {

    A <- GroupsVec[[i]]
    B <- Beta[,Pheno$Group==A]
    Medians.List[[i]] <- rowMedians(B)
    message(i)
  }

  Medians.List <- do.call(cbind,Medians.List)

  DiffList.closest <- list()

  if(verbose) {message("Generating matrix of typewise medians")}
  Medians.List <- Medians.List[complete.cases(Medians.List),]
  for(i in 1:ncol(Medians.List)) {

    M <- Medians.List
    M1 <- M[,i]
    M2 <- M[,-i]
    df <- data.frame(MinDiff= M1 - rowMin(M2),MaxDiff= M1 - rowMax(M2))
    DiffList.closest[[i]] <- df
    message(i)

  }


  DiffList2 <- lapply(DiffList.closest, function(x) transform(x,ID=rownames(x))%>%filter(MaxDiff > nnGap | MinDiff < -nnGap))

  lapply(DiffList2, dim)

  ## Intersect then

  Concise.set <- list()


  for(i in 1:length(GroupsVec)) {

    TopTab <- Filtered.tt[[i]]
    Ref <- DiffList2[[i]]$ID

    Concise.set[[i]] <- TopTab%>%filter(ID %in% Ref)
    message(i)

  }

  if(verbose) { message("Beginning final step of feature selection")}

   Concise.set2 <- lapply(Concise.set, function(x) x%>%.[order(.$t),])

   #This block of code basically sets feature upper bound - we would like symmetric feature numbers for cell types

   Threshold.val <- lapply(Concise.set2, nrow)%>%unlist(.)%>%min(.)
   if(feature.number > Threshold.val) { feature.number <- Threshold.val}
   Halfway.point <- feature.number/2
   message(feature.number)
   Sig <- lapply(Concise.set2,function(x) x[c(1:Halfway.point,(nrow(x)-Halfway.point):nrow(x)),])
   SigIDs <- do.call(rbind,Sig)%>%.$ID%>%unique(.)

   if(verbose){message("Signature Features Identified")}


   #Finally, this module selects and summarises the feature set for a run


   Beta2 <- Beta[rownames(Beta) %in% SigIDs,]
   Beta2 <- 100*Beta2
   Beta2 <- t(Beta2)
   Beta2 <- data.frame(Beta2)
   Beta2$ID <- Pheno$Group
   D2.summary <- Beta2%>%group_by(ID)%>%summarise_each(funs(mean))

   colnames(D2.summary)[[1]] <- "NAME"
   D2.summary <- t(D2.summary)
   D2.summary2 <- D2.summary[2:nrow(D2.summary),]
   colnames(D2.summary2) <- as.character(D2.summary[1,])
   D3.summary <- cbind(data.frame(NAME=rownames(D2.summary2), data.frame(D2.summary2)))
   names(Concise.set2) <- GroupsVec
   return(list(TopTables = Concise.set2, SignatureMeans = D3.summary))
   if(verbose){message("Completed")}

}

#' @export
Prep.CancerType <- function(Beta, Probes, fname) {

  Common <- intersect(Probes, rownames(Beta))
  Beta <- Beta[match(Common,rownames(Beta)),]
  NAME <- data.frame(NAME = Common)
  Beta <- cbind(NAME, data.frame(Beta))
  write.table(Beta, file = paste0(fname,".txt"),sep = "\t", row.names = FALSE, quote = F )
}



