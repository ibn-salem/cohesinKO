#!/usr/bin/Rscript
#==============================================================================================================================================
#
#     Trying out DGCA package for the cohesin KO expression data (for cherrypicking!? :) ) 
#						
#							            cohesinKO_DGCA.R             
#
#                 https://github.com/andymckenzie/DGCA
#==============================================================================================================================================

library(DGCA, quietly = TRUE)        # Load package
library(matrixStats, quietly = TRUE) # Notably, most of the filtering methods (except for rowMeans) require the matrixStats library, so we load it                                      #         here prior to filtering the genes.
library(tidyverse, quietly = TRUE)                   # Load ggplot, dplyr, etc
library(gplots, quietly = TRUE)                      # In order to visualize the differential correlation structure, DGCA offers a heatmap function. To turn it on, set the heatmapPlot argument to TRUE within ddcorAll. This requires the function heatmap.2 from the library gplots. Note that you can specify additional plotting options to the heatmap.2 at the end of the ddcorAll function.


# Load expression data for macrophages and tidy this up a bit
load("~/Dropbox/PostDoc/Rad21KOGenePairs/Data/Macrophages/RNAseqWTRad21KOMacrophages.RData")
ReadCounts.macrophages <- dplyr::select(ReadCounts.macrophages, -ensembl_gene_id)
ReadCounts.macrophages <- ReadCounts.macrophages[!duplicated( ReadCounts.macrophages[,"external_gene_id"]),]
rownames(ReadCounts.macrophages) = ReadCounts.macrophages$external_gene_id
ReadCounts.macrophages <- dplyr::select(ReadCounts.macrophages, -external_gene_id)

# Rearrange columns so that you have first the Rad21KO and then the WT
ReadCounts.macrophages <- ReadCounts.macrophages[,c(16:30, 1:15)]

# Convert to numeric; first save the rownames information
genes.names <- rownames(ReadCounts.macrophages)
ReadCounts.macrophages <- apply(ReadCounts.macrophages,2,as.numeric)
ReadCounts.macrophages <- as.data.frame(ReadCounts.macrophages)
rownames(ReadCounts.macrophages) = genes.names

# Prepare design matrix according to tutorial
n_rad21ko_samples = 15; n_wt_samples = 15
cell_type = c(rep("rad21ko", n_rad21ko_samples), rep("wt", n_wt_samples))
design_mat = model.matrix(~ cell_type + 0)
colnames(design_mat) = c("rad21ko", "wt")

# Filtering the input genes
# It is often the case that the genes with the lowest average expression levels, or the lowest variance in expression levels, are less likely to #have interesting and/or biologically relevant differences in correlation between conditions. Therefore, it is sometimes desirable to filter the #input expression matrix to remove these genes, both so that there are fewer spurious differential correlations in the resulting table, and so #that the p-value adjustment for multiple comparisons (if any) will not be affected by these genes with a lower pre-test probability of differenti#al correlation.

# In order to filter the input data in this way, DGCA offers a convenient function called filterGenes. Although the Darmanis et al. RNA data has #already been filtered to only include genes with expression levels in >= the 95th percentile in each of oligodendrocytes and neurons, we will #show how to filter this data for the purposes of demonstration.

#The first way that the input data can be filtered is by removing genes with a low average – or more precisely, a low measure of central tendency. #The default central tendency measure is the median, and the default percentile is the 30th percentile, although these can both be adjusted.

# Notably, most of the filtering methods (except for rowMeans) require the matrixStats library, so we load it here prior to filtering the genes.

nrow(ReadCounts.macrophages) #nrows before filtering
ReadCounts.macrophages.filtered = filterGenes(ReadCounts.macrophages, 
  filterTypes = "central", filterCentralType = "median", filterCentralPercentile = 0.3)
nrow(ReadCounts.macrophages.filtered)

# Finding the correlations and the significance of the correlations in each condition
cor_res = getCors(inputMat = ReadCounts.macrophages.filtered, design = design_mat)

# Finding the differential correlations between conditions
dcPairs_res = pairwiseDCor(cor_res, compare = c("rad21ko", "wt"))

# Extracting the top differentially correlated pairs from the dcPairs object
dd_pairs = dcTopPairs(dcPairs_res, nPairs = 100, classify = TRUE)

# Spearman (rank-based) differential correlation analysis
# Often times it is desirable to use a non-parametric measure for correlation because of the distribution of your input data. In this case, you #can use rank-based correlation in DGCA by setting corrType = “spearman”. A downside to using this is that it is slower to calculate the correlation matrices.

# ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat,
#  compare = c("oligodendrocyte", "neuron"),
#  adjust = "none", heatmapPlot = FALSE, nPerm = 0, corrType = "spearman", nPairs = 100)
# head(ddcor_res)

# Dissecting the differential correlation structure
# Visualization via heatmap of correlations in each condition

# Since heatmaps are best visualized with a small number of genes, we first filter the genes down to a more manageable number by selecting the    genes with the highest medians and coefficients of variation.
# Note that if you want to alter the text size in the heatmap, you can set the cexRow and cexCol arguments, as indicated above.

# macrophages_top =  filterGenes(ReadCounts.macrophages.filtered, 
#  filterTypes = c("central", "dispersion"), filterCentralPercentile = 0.75, 
#  filterDispersionPercentile = 0.75)

# ddcor_res = ddcorAll(inputMat = macrophages_top, design = design_mat,
#  compare = c("rad21ko", "wt"),
#  adjust = "none", heatmapPlot = TRUE, nPerm = 0, nPairs = "all")




#=================================================================================================================================================
# Examples of paralog genes that are in the same contact domain and show interesting, opposing, correlation patterns between the WT and the cohesinKO 
#           
#
#=================================================================================================================================================

# 1. Ms4a4d::Ms4a4b (Distance between TSSs of genes 93 kb)
plotCors(inputMat = ReadCounts.macrophages.filtered, design = design_mat,
  +     compare = c("rad21ko", "wt"), geneA = "Ms4a4d", geneB = "Ms4a4b")

# 2. Serpinf1::Serpinf2 (Distance between TSSs of genes 17 kb)
plotCors(inputMat = ReadCounts.macrophages.filtered, design = design_mat,
  +     compare = c("rad21ko", "wt"), geneA = "Serpinf1", geneB = "Serpinf2")

# 3. Ppp1r9b::Samd14 (Distance between TSSs of genes 18 kb)
plotCors(inputMat = ReadCounts.macrophages.filtered, design = design_mat,
  +     compare = c("rad21ko", "wt"), geneA = "Ppp1r9b", geneB = "Samd14")

# 4. Ms4a6b::Ms4a6c (Distance between TSSs of genes 47 kb) from almost perfect correlation to negative!
plotCors(inputMat = ReadCounts.macrophages.filtered, design = design_mat,
  +     compare = c("rad21ko", "wt"), geneA = "Ms4a6b", geneB = "Ms4a6c")

# 5. Rab7l1::5430435G22Rik (Distance between TSSs of genes ±178 kb)
plotCors(inputMat = ReadCounts.macrophages.filtered, design = design_mat,compare = c("rad21ko", "wt"), geneA = "Rab7l1", geneB = "5430435G22Rik")


#=================================================================================================================================================
# Examples of paralog genes that are in the same contact domain and are largely unaffected between the WT and the cohesinKO 
#           
#
#=================================================================================================================================================

# Slc2a6::Slc2a8
plotCors(inputMat = ReadCounts.macrophages.filtered, design = design_mat,compare = c("rad21ko", "wt"), geneA = "Slc2a6", geneB = "Slc2a8")

# Gstm7::Gstm2
plotCors(inputMat = ReadCounts.macrophages.filtered, design = design_mat,compare = c("rad21ko", "wt"), geneA = "Gstm7", geneB = "Gstm2")

# Cxcl5::Cxcl3
plotCors(inputMat = ReadCounts.macrophages.filtered, design = design_mat,compare = c("rad21ko", "wt"), geneA = "Cxcl5", geneB = "Cxcl3")


