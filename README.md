# cohesinKO

Run scripts in the following order:

1. cohesinKO.data.R
2. cohesinKO.main.R
3.
  - cohesinKO.analysis.R
  - cohesinKO.analysis_fold_change.R
  - cohesinKO.DEgenes_per_TAD.R

## Main data sets

### tidyGeneDE
contains differntial gene expression per gene for all conditions

### tidyPairsTadDE
contains gene pairs wit TAD annotation for differnt TAD groups

### tidyPairsDE
contains gene pairs with expression information

### tidyPairsTadDE
contains gene pairs with both TAD and expression information