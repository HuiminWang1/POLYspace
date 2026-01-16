# POLYspace
## Revealing multicellular neighborhoods in domain-specific contexts

<p align="center">
<img src="https://github.com/HuiminWang1/POLYspace/blob/master/POLYspace.png" width="800" />
</p>

Identifying cellular neighborhoods is essential for understanding cell–cell interactions in a spatial context. However, existing approaches often overlook the complexity of multicellular interactions and the organization of spatial domains. We present POLYspace, a general and efficient framework for discovering and analyzing cellular neighborhoods of arbitrary topology accounting for spatial domains. POLYspace formulates neighborhood identification as a subgraph searching problem and leverages C3G, a fast graph canonization algorithm we developed, to achieve scalability. Applied to one in-house dataset and three publicly available datasets spanning diverse platforms and tissues, POLYspace uncovers domain-specific cellular neighborhoods that are not captured by existing methods. These neighborhoods reveal key biological mechanisms and improve phenotype prediction.

## Broad Applicability of POLYspace
POLYspace is not limited to spatial transcriptomics or cellular neighborhood analysis.  
At its core, POLYspace implements a general optimization framework for efficient subgraph mining, together with batched graph isomorphism for the canonization of colored graphs, which naturally extends to graphlet and motif analysis.
This generality makes POLYspace broadly applicable to large-scale graph mining problems across diverse domains, including but not limited to neuroscience, social science, and chemistry.  

## Installation
You can install the released version of POLYspace from Github with the following code:

```r
# install.packages("remotes")
remotes::install_github("HuiminWang1/POLYspace")

# load the package
library(POLYspace)
```

## Dependencies 
- R version ≥ 4.1.0
- Required R packages:  
  combinat, data.table, igraph, limma, methods, parallel, partitions, RCDT, reshape2, stats, utils


## Issues
All feedback, bug reports, and suggestions are welcome.
When reporting an issue, please include a detailed, reproducible example and the output of sessionInfo() in R.


## Entering POLYspace
Details in [Tutorial](https://HuiminWang1.github.io/POLYspace/)



