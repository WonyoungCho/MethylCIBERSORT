#' @export
FeatureSelect.V4 <- function( CellLines.matrix = NULL, Heatmap = FALSE, export = TRUE, sigName, Stroma.matrix, Phenotype.stroma, deltaBeta, FDR , MaxDMRs = 1000) {

  #This function uses pairwise limma
  #This may be rather slow but whatever - it shoul do better on the basis of the nearest cell types

  require(dplyr)
  require(caret)
  require(glmnet)
  require(foreach)
  require(NMF)
  require(doParallel)
  require(matrixStats)
  require(limma)

  if(!is.null(ncol(CellLines.matrix))) {  Pheno1 <- c(rep("Cancer",ncol(CellLines.matrix)))
  Pheno2 <- c(as.character(Pheno1), as.character(Phenotype.stroma)) } else { Pheno2 <- as.character(Phenotype.stroma)}
  
  Mat2 <- Stroma.matrix
  if (!is.null(CellLines.matrix)) {
     Mat2 <- cbind(CellLines.matrix, Stroma.matrix)
  }
  message("Setting up for pairwise feature selection")


  #Work out the linear model fit here

  ContrastMatrix <- design.pairs(levels(factor(Pheno2)))
  Des <- model.matrix(~0 + Pheno2)
  colnames(Des) <- rownames(ContrastMatrix)
  Fit <- lmFit(Mat2, Des)%>%
    contrasts.fit(., ContrastMatrix)%>%
    eBayes(.)


  FitList <- list()
  for(i in 1:ncol(ContrastMatrix)) {

    FitList[[i]] <- topTable(Fit, coef = i, number = nrow(Mat2))%>%
      mutate(ID = rownames(.))%>%
      filter(adj.P.Val < FDR)

    message(paste0(i, " done"))


  }



  #Here goes the pairwise dB thing

  Transformed <- data.frame(t(Mat2))
  Split <- split(Transformed, Pheno2)
  Split <- lapply(Split, function(x) colMedians(data.matrix(x)))
  Split <- do.call(cbind, Split)
  rownames(Split) <- rownames(Mat2)

  #Then I need to actually annotate each one of these comparison topTables

  dbList <- list()
  message("Getting Delta Beta estimates")
  for(i in 1:ncol(ContrastMatrix)) {

  dB <- with(data.frame(Split), eval(parse(text = colnames(ContrastMatrix)[[i]])))
  dB <- data.frame(dB = dB, ID = rownames(Split))
  dbList[[i]] <- dB
  message(paste0 (i, " done"))
  }


  #Filter by thresholds


  dbList <- lapply(dbList, function(x) filter(x, abs(dB) > deltaBeta))
for(i in 1:length(FitList)) {

  A1 <- FitList[[i]]
  A1 <- filter(A1 , ID %in% dbList[[i]]$ID)
  A1 <- A1%>%.[rev(order(.$t)),]
  if(nrow(A1) > MaxDMRs) { A1 <-  A1[1:MaxDMRs,]                   }
  FitList[[i]] <- A1
}

  Nonzeros <- lapply(FitList, function(x) dplyr::select(x,ID))
  Nonzeros <- do.call(rbind, Nonzeros)
  Nonzeros <- filter(Nonzeros, !duplicated(ID))

  Mat3 <- Mat2[rownames(Mat2) %in% Nonzeros$ID,]

  if(Heatmap) { aheatmap(Mat3, annCol = Pheno2, labCol = NA, distfun = "euclidean", annColors = list(c("black","dodgerblue3","orange","red","cyan","firebrick3","grey","brown","dodgerblue2","yellow","purple","goldenrod","darkorange")), color = colorRampPalette(c("dodgerblue3","orange","red","firebrick3"))(200))
  }

  #Print number of selected probes
 nrow(Mat3)
 Mat3 <- 100 * Mat3

  #Then , export object
  if(export) {

    DF <- as.data.frame(t(Mat3))
    DF <- split(DF, factor(Pheno2))
    Collapsed <- lapply(DF, function(x) colMedians(data.matrix(x)))
    Collapsed <- data.frame(do.call(cbind, Collapsed))
    Collapsed <- cbind(data.frame(NAME = rownames(Mat3), stringsAsFactors = F),
                       Collapsed)

    fN <- paste0(sigName, "_Signature.txt")
    write.table(Collapsed, file = fN, sep = "\t", row.names = FALSE, quote = FALSE )

  }

  return(list(SignatureMatrix = Mat3))

}



#This function creates the pairs for the pairwise matrices
design.pairs <-
  function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    design
  }






