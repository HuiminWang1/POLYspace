# About
We developed POLYspace, a scalable and broadly applicable method for identifying and analyzing cellular neighborhoods in spatial sequencing data. POLYspace efficiently detects neighborhoods of arbitrary topology while explicitly accounting for spatial domains. It can be applied to a single sample to identify enriched or depleted neighborhoods, or across multiple samples to uncover neighborhoods associated with specific phenotypes. We applied POLYspace to four datasets, including one generated in-house, and benchmarked it against existing neighborhood inference methods. POLYspace revealed novel and biologically meaningful neighborhoods missed by other approaches, and leads to improved prediction of sample phenotypes. Together, these results demonstrate that POLYspace provides a powerful framework for disentangling complex multicellular interactions in spatial sequencing data.


## Dependencies 
- R version ≥ 4.1.0
- Required R packages:  
  combinat, data.table, igraph, limma, methods, parallel, partitions, RCDT, reshape2, stats, utils


## Issues
All feedback, bug reports, and suggestions are welcome.
When reporting an issue, please include a detailed, reproducible example and the output of sessionInfo() in R.


## Entering POLYspace
Details in [Tutorial](https://HuiminWang1.github.io/POLYspace/)