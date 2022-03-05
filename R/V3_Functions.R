#' @export
FeatureSelect.V3 <- function(Heatmap = FALSE, nCores = 4, reps.resamp = 20, export = TRUE, sigName, Stroma.matrix, Phenotype.stroma, Unlog = FALSE, BetaScale = TRUE) {

#This function uses pairwise glmnetting
#This may be rather slow but whatever - it shoul do better on the basis of the nearest cell types


require(caret)
require(glmnet)
require(foreach)
require(NMF)
require(doParallel)
require(matrixStats)

if(!is.null(ncol(CellLines.matrix))) {  Pheno1 <- c(rep("Cancer",ncol(CellLines.matrix)))
Pheno2 <- c(as.character(Pheno1), as.character(Phenotype.stroma)) } else { Pheno2 <- as.character(Phenotype.stroma)}

Mat2 <- cbind(CellLines.matrix, Stroma.matrix)

message("Setting up for pairwise feature selection")
message("Beginning elastic net selection procedure")
message("63.2 bootstrapping to be used")

#Create foreach cluster for parallelisation

Features.CVparam<- trainControl(method = "boot632",number = reps.resamp,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)

Pairs <- data.frame(t(combn(unique(Pheno2),2)), stringsAsFactors = F)
Pairs <- filter(Pairs, !X1 == X2)

#Then I do the pairwise fitting here

FitList <- list()

for(i in 1:nrow(Pairs)) {

I1 <- Phenotype.stroma == Pairs[i,]$X1 | Phenotype.stroma == Pairs[i,]$X2
M1 <- Mat2[,I1]
P1 <- as.character(Pheno2[I1])

Model <- train(x = t(M1), y = factor(P1), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0.5,1),.lambda = seq(0,0.05,by=0.01)), metric = "Kappa")

Nonzeros <- coef(Model$finalModel, s = Model$bestTune$lambda)
Nonzeros <- as.matrix(Nonzeros)
Nonzeros <- data.frame(ID = rownames(Nonzeros), Coef = as.numeric(Nonzeros[,1]))
Nonzeros <- filter(Nonzeros, !Coef == 0)
FitList[[i]] <- Nonzeros

message(paste0("pair",i," done of ", nrow(Pairs)))


  }

#Then I want to unify the nonzeros and then move on

Nonzeros <- do.call(rbind, FitList)
Nonzeros <- filter(Nonzeros, !duplicated(ID))

Mat3 <- Mat2[rownames(Mat2) %in% Nonzeros$ID,]

if(Heatmap) { aheatmap(Mat3, annCol = Pheno2, labCol = NA, distfun = "euclidean", annColors = list(c("black","dodgerblue3","orange","red","cyan","firebrick3","grey","brown","dodgerblue2","yellow","purple","goldenrod","darkorange")), color = colorRampPalette(c("dodgerblue3","orange","red","firebrick3"))(200))
}

#Print number of selected probes

if(BetaScale) { Mat3 <- 100 * Mat3}
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



