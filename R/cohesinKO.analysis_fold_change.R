#'
#' Rad21 knock-out downstream analyis of expression fold changs 
#' 

require(RColorBrewer) # for picking nice colors
require(tidyverse)    # to format, sumarize and plot tidy data

COL_GROUP = brewer.pal(9, "Set1")[c(1,9,2,4)]   # for paralog vs. sampled genes
COL_COMB = brewer.pal(12, "Paired")[c(5,6,7,8,1,2,3,4)]
pairedCols <- brewer.pal(12, "Paired")
COL_GROUP_PAIRED <- c(pairedCols[c(6,5)], "azure4", "azure3", pairedCols[c(2,1,  10,9)])
COL_ALL_GENOTYPE = brewer.pal(12, "Paired")[c(10,9)]
COL_TAD_GENOTYPE = brewer.pal(12, "Paired")[c(8,7,  4,3)] # orange, green

ALPHA = 0.3

# define some parameters
outPrefix <- "results/v03.foldChange"


# load data:
#load("results/tidyPairsTadDE.Rdata")
load("results/singleTadSourceDF.Rdata")
tidyPairsTadDE <- singleTadSourceDF

#===============================================================================
# define subset operators for tidyPairsTadDE
#===============================================================================

# select a single expression data set
singleDEcomb <- tidyPairsTadDE$cond == "WT.LPS.2"
singleDEgenotype <- tidyPairsTadDE$Treatment == "LPS" &
  tidyPairsTadDE$Time_hrs == 2

# select a single TAD data set
singleTadSource <- tidyPairsTadDE$study == "Rao2014" &
  tidyPairsTadDE$tissue == "CH12" &
  tidyPairsTadDE$TADtype == "all"

grbTadSource <- tidyPairsTadDE$study == "Rao2014" &
  tidyPairsTadDE$tissue == "CH12" &
  tidyPairsTadDE$TADtype == "GRB"

# singleTadSourceDF <-  subset(tidyPairsTadDE, singleTadSource)
# save(singleTadSourceDF, file = "results/singleTadSourceDF.Rdata")


#===============================================================================
# Plot Fold changes against each other
#===============================================================================
#-------------------------------------------------------------------------------
# use a single Genotype
#-------------------------------------------------------------------------------
subDF <- subset(tidyPairsTadDE, singleDEcomb & singleTadSource)

# p <- ggplot(subDF, aes(x = log2FoldChange_1, y = log2FoldChange_2)) + 
#   geom_point(alpha = 0.1)

p <- ggplot(subDF, aes(x = log2FoldChange_1, y = log2FoldChange_2)) +
  stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  theme_bw() +
  theme(text = element_text(size = 15)) + 
  labs(x = "log2 fold-change gene_1", y =  "log2 fold-change gene_2")

ggsave(p, file = paste0(
  outPrefix, 
  ".all_groups.lfc_g1_vs_g2.WT_LPS_2.scatter.pdf"), w = 6, h = 4)

#-------------------------------------------------------------------------------
# use a single condition by Genotype
#-------------------------------------------------------------------------------
subDF <- subset(tidyPairsTadDE, singleDEgenotype & singleTadSource)

# p <- ggplot(subDF, aes(x = log2FoldChange_1, y = log2FoldChange_2)) + 
#   geom_point(alpha = 0.1)

p <- ggplot(subDF, aes(x = log2FoldChange_1, y = log2FoldChange_2)) +
  stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  facet_grid(Genotype ~ .) + 
  theme_bw() +
  theme(text = element_text(size = 15)) + 
  labs(x = "log2 fold-change gene_1", y =  "log2 fold-change gene_2")

ggsave(p, file = paste0(
  outPrefix, 
  ".all_groups.lfc_g1_vs_g2_by_Genotype.LPS_2.scatter.pdf"), w = 6, h = 6)


#===============================================================================
# compare fold changes by TAD and genotype using boxplots
#===============================================================================

# use a single condition
subDF <- subset(tidyPairsTadDE, singleDEgenotype & singleTadSource)

# %>% 
# filter(abs(log2FoldChange_1) >= .5 & abs(log2FoldChange_2) >= .5)

p <- ggplot(subDF, aes(x = Genotype, color = cond, y = lfc_diff)) + 
  geom_boxplot() + 
  facet_grid(. ~ Boundary) + 
  theme_bw() + 
  theme(text = element_text(size = 15), 
      axis.text.x = element_text(angle = 45, hjust = 1), 
      legend.position = "none") + 
  scale_color_manual(values = COL_ALL_GENOTYPE, guide_legend(title = "")) +
  guides(fill = guide_legend(title = "")) +
  labs(y = "Fold change difference\n abs(lfc_g1 - lfc_g2)", x = "")
p
  
ggsave(p, file = paste0(
  outPrefix, 
  ".all_groups.lfc_diff_by_genotype_Boundary_Rao2014_CH12.LPS_2.boxplot.pdf"), 
  w = 4, h = 6)

#-------------------------------------------------------------------------------
# use a single TAD data set but all expression combinations
#-------------------------------------------------------------------------------
subDF <- subset(tidyPairsTadDE, singleTadSource)

p <- ggplot(subDF, aes(x = Boundary, color = Genotype, y = lfc_diff)) + 
  geom_boxplot() + 
  facet_grid(Treatment ~ Time_hrs) + 
  theme_bw() + 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "bottom") + 
  scale_color_manual(values = COL_ALL_GENOTYPE, guide_legend(title = "")) +
  guides(fill = guide_legend(title = "")) +
  labs(y = "Fold change difference\n abs(lfc_g1 - lfc_g2)", x = "")

ggsave(p, file = paste0(
  outPrefix, 
  ".all_groups.lfc_diff_by_Genotype_Boundary_Time_Treatment.Rao2014.boxplot.pdf"), 
  w = 6, h = 6)

#===============================================================================
# compare fold changes by TAD and genotype by distance with scatter plot
#===============================================================================

subDF <- subset(tidyPairsTadDE, singleDEcomb & singleTadSource)

p <- ggplot(subDF, aes(x = dist, y = lfc_diff)) +
  stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  theme_bw() +
  theme(text = element_text(size = 15)) + 
  labs(x = "Distance between gene pairs [kb]",
       y = "Fold change difference\n abs(lfc_g1 - lfc_g2)")

ggsave(p, file = paste0(
  outPrefix, 
  ".all_groups.lfc_diff_by_distance.WT_LPS_2.scatter.pdf"), w = 6, h = 4)
#-------------------------------------------------------------------------------
# compare by genotype
#-------------------------------------------------------------------------------
subDF <- subset(tidyPairsTadDE, singleDEgenotype & singleTadSource)

p <- ggplot(subDF, aes(x = dist, y = lfc_diff)) +
  stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE) +
  scale_fill_continuous(low = "white", high = "dodgerblue4") +
  facet_grid(. ~ Genotype) +
  theme_bw() +
  theme(text = element_text(size = 15)) + 
  labs(x = "Distance between gene pairs [kb]",
       y = "Fold change difference\n abs(lfc_g1 - lfc_g2)")

ggsave(p, file = paste0(
  outPrefix, 
  ".all_groups.lfc_diff_by_distance_Genotype.LPS_2.scatter.pdf"), w = 6, h = 4)



#===============================================================================
# compare fold changes by TAD and genotype by distance
#===============================================================================

# use a single condition
subDF <- subset(tidyPairsTadDE, singleDEgenotype & singleTadSource)

p <- ggplot(subDF, aes(x = dist, y = lfc_diff, 
                       linetype = Genotype, fill = Boundary, col = Boundary)) + 
  geom_smooth(alpha = ALPHA) + 
  guides(col = guide_legend(title = ""), fill = guide_legend(title = ""),
         linetype = guide_legend(title = "")) + 
  scale_color_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + 
  scale_fill_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + 
  theme_bw() + 
  theme(legend.position = "bottom") +
  labs(
    x = "Distance between gene pairs [kb]",
    y = "Fold change difference\n abs(lfc_g1 - lfc_g2)")

ggsave(p, file = paste0(
  outPrefix, 
  ".all_groups.lfc_diff_by_distance_genotype_Boundary_Rao2014_CH12.LPS_2.smooth.pdf"), 
  w = 6, h = 6)

#-------------------------------------------------------------------------------
# use a single TAD data set but all expression conditions
#-------------------------------------------------------------------------------
subDF <- subset(tidyPairsTadDE, singleTadSource)

p <- ggplot(subDF, aes(x = dist, y = lfc_diff,
                         linetype = Genotype, fill = Boundary, col = Boundary)) +
  geom_smooth(alpha = ALPHA) +
  facet_grid(Treatment ~ Time_hrs) +
  guides(col = guide_legend(title = ""), fill = guide_legend(title = ""),
         linetype = guide_legend(title = "")) +
  scale_color_manual(values = COL_TAD_GENOTYPE[c(1,3)]) +
  scale_fill_manual(values = COL_TAD_GENOTYPE[c(1,3)]) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    x = "Distance between gene pairs [kb]",
    y = "Fold change difference\n abs(lfc_g1 - lfc_g2)")

ggsave(p, file = paste0(
  outPrefix,
  ".all_groups.lfc_diff_by_distance_genotype_Boundary_treatment_Rao2014_CH12.LPS_2.smooth.pdf"),
  w = 6, h = 6)



#===============================================================================
# Ratio of fold changes 
#===============================================================================


#===============================================================================
# Plot Fold changes against each other
#===============================================================================

#-------------------------------------------------------------------------------
# compare fold changes ratio by boxplot
#-------------------------------------------------------------------------------

# use a single condition
subDF <- subset(tidyPairsTadDE, singleDEgenotype & singleTadSource)

p <- ggplot(subDF, aes(x = Genotype, color = cond, y = (log2FoldChange_1 / log2FoldChange_2))) + 
  geom_boxplot() + 
  facet_grid(. ~ Boundary) + 
  theme_bw() + 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") + 
  scale_color_manual(values = COL_ALL_GENOTYPE, guide_legend(title = "")) +
  guides(fill = guide_legend(title = "")) +
  labs(y = "Fold change ratio \n lfc_g1  /  lfc_g2", x = "")

ggsave(p, file = paste0(
  outPrefix, 
  ".all_groups.lfc_ratio_by_genotype_Boundary_Rao2014_CH12.LPS_2.boxplot.pdf"), 
  w = 4, h = 6)


#-------------------------------------------------------------------------------
# compare fold changes log ratio by boxplot
#-------------------------------------------------------------------------------

# use a single condition
subDF <- subset(tidyPairsTadDE, singleDEgenotype & singleTadSource)

p <- ggplot(subDF, aes(x = Genotype, color = cond, y = log2(log2FoldChange_1 / log2FoldChange_2))) + 
  geom_boxplot() + 
  facet_grid(. ~ Boundary) + 
  theme_bw() + 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none") + 
  scale_color_manual(values = COL_ALL_GENOTYPE, guide_legend(title = "")) +
  guides(fill = guide_legend(title = "")) +
  labs(y = "Fold change log ratio \n log2(lfc_g1  /  lfc_g2)", x = "")

ggsave(p, file = paste0(
  outPrefix, 
  ".all_groups.lfc_logRatio_by_genotype_Boundary_Rao2014_CH12.LPS_2.boxplot.pdf"), 
  w = 4, h = 6)



#===============================================================================
# compare fold changes log ratio by TAD and genotype by distance
#===============================================================================

# use a single condition
subDF <- subset(tidyPairsTadDE, singleDEgenotype & singleTadSource)

p <- ggplot(subDF, aes(x = dist, y = log2(log2FoldChange_1 / log2FoldChange_2), 
                       linetype = Genotype, fill = Boundary, col = Boundary)) + 
  geom_smooth(alpha = ALPHA) + 
  guides(col = guide_legend(title = ""), fill = guide_legend(title = ""),
         linetype = guide_legend(title = "")) + 
  scale_color_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + 
  scale_fill_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + 
  theme_bw() + 
  theme(legend.position = "bottom") +
  labs(
    x = "Distance between gene pairs [kb]",
    y = "Fold change log ratio \n log2(lfc_g1  /  lfc_g2)")

ggsave(p, file = paste0(
  outPrefix, 
  ".all_groups.lfc_logRatio_by_distance_genotype_Boundary_Rao2014_CH12.LPS_2.smooth.pdf"), 
  w = 6, h = 6)

