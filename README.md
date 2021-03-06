## Quick link

| Function | Code |
| -------- | ---- |
| V4 | https://github.com/WonyoungCho/MethylCIBERSORT/blob/main/R/V4_Functions.R |
| toptable | https://github.com/WonyoungCho/MethylCIBERSORT/blob/main/R_ref/toptable.R |
| lmfit | https://github.com/WonyoungCho/MethylCIBERSORT/blob/main/R_ref/lmfit.R |

---

# MethylCIBERSORT
[Pan-cancer deconvolution of tumour composition using DNA methylation (2018)](https://www.nature.com/articles/s41467-018-05570-1)

MethylCIBERSORT v2.0.1 https://zenodo.org/record/1298968#.YiMc5ehBzmg

```bash
R CMD INSTALL MethylCIBERSORT_0.2.1.tar.gz
```

```R
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
beta_file <- args[1]
project <- args[2]

## Mat <- read.table(file="beta_values.tsv", sep = '\t', header = TRUE, row.names=1)
Mat <- read.table(file=beta_file, sep = '\t', header = TRUE, row.names=1)
head(Mat)

library("MethylCIBERSORT")

data("StromalMatrix_V2")

Stromal_v2 <- Stromal_v2[sort(rownames(Stromal_v2)),]

Int <- intersect(rownames(Mat), rownames(Stromal_v2))
Mat <- Mat[match(Int, rownames(Mat)),]
Stromal_v2 <- Stromal_v2[match(Int, rownames(Stromal_v2)),]

RefData <- Stromal_v2
RefPheno <- Stromal_v2.pheno

Signature <- FeatureSelect.V4(CellLines.matrix = NULL,
                              Heatmap = TRUE,
                              export = TRUE,
                              sigName = "MyReference",
                              Stroma.matrix = RefData,
                              deltaBeta = 0.2,
                              FDR = 0.01,
                              MaxDMRs = 100,
                              Phenotype.stroma = RefPheno)
```

# Cell type
- CD4+  : regulatory T (Treg) cells and conventional T helper (Th) cells.
- CD14+ : Monocyte and macrophages.
- CD19+ : B cells
- CD56+ : NK cells


Reference
- https://github.com/dannlbol/mcibersort_scripts
- https://github.com/cran/limma/tree/master/R
- 
