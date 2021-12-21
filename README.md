# GOfisher  

This is an R packge for simple GO enrichment analyses of *Arabidopsis thaliana*.  
The GOfisher package includes a sort of functions to perform Fisher's exact probability tests with FDR correction.  

## Installation
The package is simply installed via this GitHub repository using the devtools package.  
```
devtools::install_github("https://github.com/yassato/GOfisher")
```

## Usage

### Input files  
A list of AGI codes and GO terms is made from [the TAIR annotation file](https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt)
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
