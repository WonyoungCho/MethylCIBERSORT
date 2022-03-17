# Script.R

```R
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
beta_file <- args[1]
project <- args[2]

library(data.table)

Mat <- data.frame(fread(beta_file), row.names=1)
Stromal_v2 <- data.frame(fread("Stromal_v3.txt"), row.names=1)

Int <- intersect(rownames(Mat), rownames(Stromal_v2))
Mat <- Mat[match(Int, rownames(Mat)),]
Stromal_v2 <- Stromal_v2[match(Int, rownames(Stromal_v2)),]


RefData <- Stromal_v2
RefPheno <- readRDS("Stromal_v3.RData")


if(!exists("foo", mode="function")) source("V4_Functions.R")
Signature <- FeatureSelect.V4(CellLines.matrix = NULL,
                              Heatmap = FALSE,
                              export = TRUE,
                              sigName = project,
                              Stroma.matrix = RefData,
                              deltaBeta = 0.2,
                              FDR = 0.01,
                              MaxDMRs = 100,
                              Phenotype.stroma = RefPheno)
```
