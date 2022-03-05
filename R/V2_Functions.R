#' @export

Featureselect.V2 <- function(CellLines.matrix, Heatmap = FALSE, nCores = 4, reps.resamp = 20, export = TRUE, sigName, Stroma.matrix, Phenotype.stroma, Unlog = TRUE) {

  require(caret)
  require(glmnet)
  require(foreach)
  require(NMF)
  require(doParallel)
  require(matrixStats)

if(!is.null(ncol(CellLines.matrix))) {  Pheno1 <- c(rep("Cancer",ncol(CellLines.matrix)))
  Pheno2 <- c(as.character(Pheno1), as.character(Phenotype.stroma)) }

  else { Pheno2 <- as.character(Phenotype.stroma)}

  Mat2 <- cbind(CellLines.matrix, Stroma.matrix)

  message("Beginning elastic net selection procedure")
  message("63.2 bootstrapping to be used")

  #Create foreach cluster for parallelisation

Features.CVparam<- trainControl(method = "boot632",number = reps.resamp,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)

if(nCores > 1) {  registerDoParallel(makeCluster(nCores))
  message( "Parallelisation schema set up")}


Model <- train(x = t(Mat2), y = factor(Pheno2), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0.5,1),.lambda = seq(0,0.05,by=0.01)), metric = "Kappa")

message("Retrieving Nonzero Coefficients")
Nonzeros <-  coef(Model$finalModel, s = Model$bestTune$lambda)
Nonzeros <- lapply(Nonzeros, function(x) data.frame(ID = rownames(x), Coef = as.numeric(x[,1])))
Nonzeros <- lapply(Nonzeros, function(x) filter(x, !Coef == 0))
Nonzeros <- do.call(rbind, Nonzeros)

#Then I need to do the whole shebang of getting the features

Mat3 <- Mat2[rownames(Mat2) %in% Nonzeros$ID,]

if(Heatmap) { aheatmap(Mat3, annCol = Pheno2, labCol = NA, distfun = "euclidean", annColors = list(c("black","dodgerblue3","orange","red","cyan","firebrick3","grey","brown","dodgerblue2","yellow","purple","goldenrod","darkorange")), color = colorRampPalette(c("dodgerblue3","orange","red","firebrick3"))(200))
  }

#Print number of selected probes

print(nrow(Mat3))
if(Unlog) { Mat3 <- 2^Mat3}


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



