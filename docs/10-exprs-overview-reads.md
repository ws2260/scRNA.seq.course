---
output: html_document
---

## Data visualization (Reads)


```r
library(scater)
options(stringsAsFactors = FALSE)
reads <- readRDS("tung/reads.rds")
reads.qc <- reads[rowData(reads)$use, colData(reads)$use]
endog_genes <- !rowData(reads.qc)$is_feature_control
```




```r
plotPCA(
    reads[endog_genes, ],
    exprs_values = "counts",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```

<div class="figure" style="text-align: center">
<img src="10-exprs-overview-reads_files/figure-html/expr-overview-pca-before-qc-reads1-1.png" alt="PCA plot of the tung data" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-before-qc-reads1)PCA plot of the tung data</p>
</div>


```r
plotPCA(
    reads[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```

<div class="figure" style="text-align: center">
<img src="10-exprs-overview-reads_files/figure-html/expr-overview-pca-before-qc-reads2-1.png" alt="PCA plot of the tung data" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-before-qc-reads2)PCA plot of the tung data</p>
</div>


```r
plotPCA(
    reads.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```

<div class="figure" style="text-align: center">
<img src="10-exprs-overview-reads_files/figure-html/expr-overview-pca-after-qc-reads-1.png" alt="PCA plot of the tung data" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-after-qc-reads)PCA plot of the tung data</p>
</div>


```r
plotTSNE(
    reads[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 130,
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    rand_seed = 123456
)
```

<div class="figure" style="text-align: center">
<img src="10-exprs-overview-reads_files/figure-html/expr-overview-tsne-before-qc-reads-1.png" alt="tSNE map of the tung data" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-before-qc-reads)tSNE map of the tung data</p>
</div>


```r
plotTSNE(
    reads.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 130,
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    rand_seed = 123456
)
```

<div class="figure" style="text-align: center">
<img src="10-exprs-overview-reads_files/figure-html/expr-overview-tsne-after-qc-reads-1.png" alt="tSNE map of the tung data" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-after-qc-reads)tSNE map of the tung data</p>
</div>

<div class="figure" style="text-align: center">
<img src="10-exprs-overview-reads_files/figure-html/expr-overview-tsne-after-qc-exercise2-1-1.png" alt="tSNE map of the tung data (perplexity = 10)" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-after-qc-exercise2-1)tSNE map of the tung data (perplexity = 10)</p>
</div>

<div class="figure" style="text-align: center">
<img src="10-exprs-overview-reads_files/figure-html/expr-overview-tsne-after-qc-exercise2-2-1.png" alt="tSNE map of the tung data (perplexity = 200)" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-after-qc-exercise2-2)tSNE map of the tung data (perplexity = 200)</p>
</div>


```r
sessionInfo()
```

```
## R version 3.4.2 (2017-09-28)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 9 (stretch)
## 
## Matrix products: default
## BLAS/LAPACK: /usr/lib/libopenblasp-r0.2.19.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] knitr_1.17                 scater_1.6.0              
##  [3] SingleCellExperiment_1.0.0 SummarizedExperiment_1.8.0
##  [5] DelayedArray_0.4.1         matrixStats_0.52.2        
##  [7] GenomicRanges_1.30.0       GenomeInfoDb_1.14.0       
##  [9] IRanges_2.12.0             S4Vectors_0.16.0          
## [11] ggplot2_2.2.1              Biobase_2.38.0            
## [13] BiocGenerics_0.24.0       
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.4.0           edgeR_3.20.1           
##  [3] bit64_0.9-7             viridisLite_0.2.0      
##  [5] shiny_1.0.5             assertthat_0.2.0       
##  [7] highr_0.6               blob_1.1.0             
##  [9] GenomeInfoDbData_0.99.1 vipor_0.4.5            
## [11] yaml_2.1.14             progress_1.1.2         
## [13] RSQLite_2.0             backports_1.1.1        
## [15] lattice_0.20-34         glue_1.2.0             
## [17] limma_3.34.0            digest_0.6.12          
## [19] XVector_0.18.0          colorspace_1.3-2       
## [21] cowplot_0.8.0           htmltools_0.3.6        
## [23] httpuv_1.3.5            Matrix_1.2-7.1         
## [25] plyr_1.8.4              XML_3.98-1.9           
## [27] pkgconfig_2.0.1         biomaRt_2.34.0         
## [29] bookdown_0.5            zlibbioc_1.24.0        
## [31] xtable_1.8-2            scales_0.5.0           
## [33] Rtsne_0.13              tibble_1.3.4           
## [35] lazyeval_0.2.1          magrittr_1.5           
## [37] mime_0.5                memoise_1.1.0          
## [39] evaluate_0.10.1         beeswarm_0.2.3         
## [41] shinydashboard_0.6.1    tools_3.4.2            
## [43] data.table_1.10.4-3     prettyunits_1.0.2      
## [45] stringr_1.2.0           munsell_0.4.3          
## [47] locfit_1.5-9.1          AnnotationDbi_1.40.0   
## [49] bindrcpp_0.2            compiler_3.4.2         
## [51] rlang_0.1.4             rhdf5_2.22.0           
## [53] grid_3.4.2              RCurl_1.95-4.8         
## [55] tximport_1.6.0          rjson_0.2.15           
## [57] labeling_0.3            bitops_1.0-6           
## [59] rmarkdown_1.7           gtable_0.2.0           
## [61] DBI_0.7                 reshape2_1.4.2         
## [63] R6_2.2.2                gridExtra_2.3          
## [65] dplyr_0.7.4             bit_1.1-12             
## [67] bindr_0.1               rprojroot_1.2          
## [69] stringi_1.1.5           ggbeeswarm_0.6.0       
## [71] Rcpp_0.12.13
```
