# GOfisher  

This is an R packge for simple GO enrichment analyses on *Arabidopsis thaliana* genes.  
The GOfisher package includes a sort of functions to perform Fisher's exact probability tests with FDR correction.  
[![DOI](https://zenodo.org/badge/440548100.svg)](https://zenodo.org/badge/latestdoi/440548100)  

## Installation
The package can be installed via this GitHub repository using the devtools package.  
```
devtools::install_github("https://github.com/yassato/GOfisher")
```

## Usage

### Input files  
A list of AGI codes and GO terms is made from [the TAIR annotation files](https://doi.org/10.5281/zenodo.7159104).  
```
fn <- "./PATH/ATH_GO_GOSLIM.txt"
ulg <- ng.GOprep_TAIR(fn) # it takes ca. 3 hrs
save(ulg,file="ulg", version=2)
```

### Run Fisher tests
Fisher tests are run multiple times using ng.mft(), and then corrected by the false discovery rate using ng.prepGOtestOutTable().
```
data(ulg)
data(gl)
fisher.res <- ng.mft(ulg, gl)
GO.list <- ng.prepGOtestOutTable(fisher.res[fisher.res[,"xtt"]>1,], alpha=0.05)
GO.list[order(GO.list[,1]),]
```

### Visualization
The GOfisher package has the utility function, named GOfisher2REVIGO(), to convert its output as an input of the rrvgo package.
```
library(rrvgo)
library(org.At.tair.db)
res <- GOfisher2REVIGO(fisher.res[fisher.res[,"xtt"]>1,],gl,ulg)
simMatrix <- calculateSimMatrix(res$ID, orgdb="org.At.tair.db", ont="BP", method="Rel")
scores <- setNames(-log10(res$qvalue), res$ID)
reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.At.tair.db")

treemapPlot(reducedTerms)
```

### Tips  
To search genes with a given GO ID  
```
data(ulg)
data(gl)
genes <- ng.SearchGOTerms(gl, "GO:0006952", ulg)
```

Then to check gene descriptions  
```
data(des)
des[genes,]
```

## References
The original source codes are developed by [A.J. Nagano Lab](https://github.com/naganolab/). The same methods were used by Kamitani et al. (2019); Sato et al. (2019); and others.  
- Kamitani M, Kashima M, Tezuka A, Nagano AJ. (2019) Lasy-Seq: A High-Throughput Library Preparation Method for RNA-Seq and Its Application in the Analysis of Plant Responses to Fluctuating Temperatures. Scientific Reports 9:7091. https://doi.org/10.1038/s41598-019-43600-0.  
- Sato Y, Tezuka A, Kashima M, Deguchi A, Shimizu-Inatsugi R, Yamazaki M, Shimizu KK, Nagano AJ. (2019) Transcriptional Variation in Glucosinolate Biosynthetic Genes and Inducible Responses to Aphid Herbivory on Field-Grown *Arabidopsis Thaliana*. Frontiers in Genetics 10:787. https://doi.org/10.3389/fgene.2019.00787.


