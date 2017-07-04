# cohesinKO

We hope this to be a package containing functions for analyzing cohesin KO data

## Script

Run scripts in the following order:

1. [cohesinKO.data.R](R/cohesinKO.data.R)
2. [cohesinKO.main.R](R/cohesinKO.main.R)
3. one or all of the following:
  - [cohesinKO.analysis.R](R/cohesinKO.analysis.R)
  - [cohesinKO.analysis_fold_change.R](R/cohesinKO.analysis_fold_change.R)
  - [cohesinKO.DEgenes_per_TAD.R](R/cohesinKO.DEgenes_per_TAD.R)

## Main data sets

the main dataset which are created and used by the scripts are.


- **tidyGeneDE** contains differential gene expression per gene for all conditions

- **tidyPairsTadDE** contains gene pairs with TAD annotation for different TAD groups

- **tidyPairsDE** contains gene pairs with expression information

- **tidyPairsTadDE** contains gene pairs with both TAD and expression information
