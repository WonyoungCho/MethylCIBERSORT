# MethylCIBERSORT
[Pan-cancer deconvolution of tumour composition using DNA methylation](https://www.nature.com/articles/s41467-018-05570-1)

MethylCIBERSORT v2.0.1 https://zenodo.org/record/1298968#.YiMc5ehBzmg

```bash
R CMD INSTALL MethylCIBERSORT_0.2.1.tar.gz
```

```R
library("MethylCIBERSORT")

Mat <- read.table(file = 'beta_values.tsv', sep = '\t', header = TRUE, row.names = 1)

data("StromalMatrix_V2")

Stromal_v2 <- Stromal_v2[sort(rownames(Stromal_v2)),]

Int <- intersect(rownames(Mat), rownames(Stromal_v2))
Mat <- Mat[match(Int, rownames(Mat)),]
Stromal_v2 <- Stromal_v2[match(Int, rownames(Stromal_v2)),]

RefData <- Stromal_v2
RefPheno <- Stromal_v2.pheno

Signature <- FeatureSelect.V4(CellLines.matrix = NULL,
                              Heatmap = FALSE,
                              export = TRUE,
                              sigName = "MyReference",
                              Stroma.matrix = RefData,
                              deltaBeta = 0.2,
                              FDR = 0.01,
                              MaxDMRs = 100,
                              Phenotype.stroma = RefPheno)
```


Reference
- https://github.com/dannlbol/mcibersort_scripts
