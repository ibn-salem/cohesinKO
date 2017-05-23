#'
#' Analyse the number of responsive genes (differntial expressed by treatment) 
#' per TAD.
#' 

require(RColorBrewer) # for picking nice colors
require(tidyverse)    # to format, sumarize and plot tidy data

COL_COMB = brewer.pal(12, "Paired")[c(5,6,7,8,1,2,3,4)]
COL_ALL_GENOTYPE = brewer.pal(12, "Paired")[c(10,9)]
COL_TAD_GENOTYPE = brewer.pal(12, "Paired")[c(8,7,  4,3)] # orange, green
COL_RAND = brewer.pal(8, "Accent")[c(7, 8)] # bown, gray

ALPHA = 0.3

# define some parameters
outPrefix <- "results/v03.DEgenes_per_TAD"


# load data:
load("results/TAD2geneDF.Rdata")
load("results/tidyGeneDE.Rdata")
load("results/genesGR.Rdata")

#===============================================================================
# Add DE to TAD2genesDF
#===============================================================================

tidyTAD2gene <- TAD2geneDF %>% 
  left_join(tidyGeneDE, by = c("ENSG" = "EnsemblID"))

#-------------------------------------------------------------------------------
# vectors for subsetting the data set
#-------------------------------------------------------------------------------
# select a single expression data set
singleDEcomb <- tidyTAD2gene$cond == "WT.LPS.8"

singleDEcondition <- tidyTAD2gene$Treatment == "LPS" &
  tidyTAD2gene$Time_hrs == 8


# select a single TAD origin
singleTadorigen <- tidyTAD2gene$study == "Rao2014" &
  tidyTAD2gene$tissue == "CH12" 

# select a single TAD data set
singleTadSource <- tidyTAD2gene$study == "Rao2014" &
  tidyTAD2gene$tissue == "CH12" &
  tidyTAD2gene$GRB == "all"

#===============================================================================
# DE genes by condition
#===============================================================================

# get subset of genes that are contained in tidyGeneDE
subTidyGeneDE <- tidyGeneDE %>% 
  filter(EnsemblID %in% genesGR$gene_id)


countGeneDF <- subTidyGeneDE %>% 
  count(Genotype, Treatment, Time_hrs, DE) 

p <- ggplot(countGeneDF, aes(x = Genotype, y = n, fill = DE)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(1), vjust = 1.5) + 
  facet_grid(Treatment ~ Time_hrs) +
  theme_bw() +
  theme(text = element_text(size = 15)) + 
  labs(y = "Genes")

ggsave(p, file = paste0(outPrefix, 
                        ".genes_by_DE_genotype_conditions.pdf"), w = 6, h = 6)

#===============================================================================
# Analyse number of genes per TAD
#===============================================================================

subDF <- subset(tidyTAD2gene, singleDEcomb & singleTadorigen)

countDF <- subDF %>%
  group_by(GRB, TADid) %>% 
  summarize(nGenes = n()) %>% 
  group_by(GRB, nGenes) %>% 
  summarize(n = n()) %>% 
  mutate(percent = n / sum(n) * 100)

p <- ggplot(countDF, aes(x = nGenes, y = percent)) +
  geom_bar(stat = "identity") +
  facet_grid(GRB ~ .)  + 
  theme_bw() +
  theme(text = element_text(size = 15)) + 
  labs(x = "Genes per TAD", y =  "%")

ggsave(p, file = paste0(outPrefix, 
  ".genes_per_TAD_by_GRB.pdf"), w = 6, h = 6)


#-------------------------------------------------------------------------------
# Analyse DE gens per TAD
#-------------------------------------------------------------------------------
subDF <- subset(tidyTAD2gene, singleDEcomb & singleTadorigen)

countDF <- subDF %>%
  group_by(TADid, GRB, DE) %>% 
  summarize(nGenes = n()) %>% 
  mutate(nGenes = factor(
    ifelse(nGenes < 10, nGenes, ">= 10"),
    c(parse_character(1:9), ">= 10")
  )) %>% 
  group_by(GRB, DE, nGenes) %>% 
  summarize(n = n()) %>% 
  mutate(percent = n / sum(n) * 100)

p <- ggplot(countDF, aes(x = nGenes, y = percent, fill = DE)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(GRB ~ DE)  + 
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "none") + 
  labs(x = "Genes per TAD", y =  "%")


ggsave(p, file = paste0(outPrefix, 
                        ".DEgenes_per_TAD_by_GRB.WT_LPS_8.pdf"), w = 6, h = 6)

#-------------------------------------------------------------------------------
# Analyse DE gens per TAD by Genotype
#-------------------------------------------------------------------------------
subDF <- subset(tidyTAD2gene, singleDEcondition & singleTadorigen)

countDF <- subDF %>%
  group_by(TADid, GRB, Genotype, DE) %>% 
  summarize(nGenes = n()) %>% 
  mutate(nGenes = factor(
    ifelse(nGenes < 10, nGenes, "10+"),
    c(parse_character(1:9), "10+")
  )) %>% 
  mutate(nGenesBinary = factor(
    ifelse(nGenes == 1, "1", ">1")
  ))

plotDF <- countDF %>% 
  group_by(GRB, Genotype, DE, nGenes) %>% 
  summarize(n = n()) %>% 
  mutate(percent = n / sum(n) * 100)

# calculate fisher test on each group
pValDF <- countDF %>% 
  group_by(GRB, DE) %>% 
  summarize(
    p = fisher.test(nGenesBinary, Genotype)$p.value
  )

p <- ggplot(plotDF, aes(x = nGenes, y = percent)) +
  geom_bar(aes(fill = Genotype), stat = "identity", position = "dodge") +
  geom_text(aes(x = 5, y = 70, label = paste0("p = ", signif(p, 3))), data = pValDF) +
  facet_grid(GRB ~ DE)  + 
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "bottom") + 
  labs(x = "Genes per TAD", y =  "% of TADs")

ggsave(p, file = paste0(outPrefix, 
                        ".DEgenes_per_TAD_by_GRB_and_Genotype.Rao2014_CH12.LPS_8.pdf"), w = 6, h = 6)

#-------------------------------------------------------------------------------
# Analyse DE gens per TAD compared to randomized
#-------------------------------------------------------------------------------


subDF <- subset(tidyTAD2gene, singleDEcomb & singleTadorigen)

# check number of DE genes per 
# subDF %>% count(GRB, DE)

# add randmized data column per GRB group
realRandDF <- subDF %>% 
  group_by(GRB) %>% 
  do(mutate(., randDE = sample(DE, replace = TRUE))) %>% 
  ungroup() %>% 
  gather(DE, randDE, key = "expGroup", value = "DE") %>% 
  mutate(expGroup = ifelse(expGroup == "DE", "actual", "permuted"))

# subDF %>% count(GRB, randDE)

countDF <- realRandDF %>%
  group_by(TADid, GRB, expGroup, DE) %>% 
  summarize(nGenes = n()) %>% 
  mutate(nGenes = factor(
    ifelse(nGenes <= 9, nGenes, "10+"),
    c(parse_character(1:9), "10+")
  )) %>% 
  mutate(nGenesBinary = factor(
    ifelse(nGenes == 1, "1", ">1")
  ))

# calculate fisher test on each group
pValDF <- countDF %>% 
  group_by(GRB, DE) %>% 
  summarize(
    p = fisher.test(nGenesBinary, expGroup)$p.value
  )

plotDF <- countDF %>%   
  group_by(GRB, expGroup, DE, nGenes) %>% 
  summarize(n = n()) %>% 
  mutate(percent = n / sum(n) * 100)

p <- ggplot(plotDF, aes(x = nGenes, y = percent)) +
  geom_bar(aes(fill = expGroup), stat = "identity", position = "dodge") +
  # geom_text(aes(label = n), angle = 90, hjust = 1.6, position = position_dodge(.9), size = 3) +
  geom_text(aes(x = 5, y = 70, label = paste0("p = ", signif(p, 3))), data = pValDF) +
  facet_grid(GRB ~ DE)  + 
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "bottom") + 
  scale_fill_manual(values = COL_RAND) +
  labs(x = "Genes per TAD", y =  "% of TADs")

ggsave(p, file = paste0(outPrefix, 
                        ".DEgenes_per_TAD_by_GRB_and_permutation.Rao2014_CH12.LPS_8.pdf"), w = 6, h = 6)


#-------------------------------------------------------------------------------
# Analyse DE gens per TAD compared to randomized fro all TAD data sets
#-------------------------------------------------------------------------------
subDF <- subset(tidyTAD2gene, singleDEcomb)

# check number of DE genes per 
# subDF %>% count(GRB, DE)

# add randmized data column per GRB group
realRandDF <- subDF %>% 
  group_by(study, tissue, GRB) %>% 
  do(mutate(., randDE = sample(DE, replace = TRUE))) %>% 
  ungroup() %>% 
  gather(DE, randDE, key = "expGroup", value = "DE") %>% 
  mutate(expGroup = ifelse(expGroup == "DE", "actual", "permuted"))

# subDF %>% count(GRB, randDE)

countDF <- realRandDF %>%
  group_by(TADid, study, tissue, GRB, expGroup, DE) %>% 
  summarize(nGenesInt = n()) %>% 
  mutate(nGenes = factor(
    ifelse(nGenesInt <= 9, nGenesInt, "10+"),
    c(parse_character(1:9), "10+")
  )) %>% 
  mutate(nGenesBinary = factor(
    ifelse(nGenesInt == 1, "1", ">1")
  ))

# calculate fisher test on each group
pValDF <- countDF %>% 
  group_by(study, tissue, GRB, DE) %>% 
  summarize(
    p = fisher.test(nGenesBinary, expGroup)$p.value
  )

plotDF <- countDF %>%   
  group_by(study, tissue, GRB, expGroup, DE, nGenesBinary) %>% 
  summarize(n = n()) %>% 
  mutate(percent = n / sum(n) * 100)

p <- ggplot(plotDF, aes(x = nGenesBinary, y = percent)) +
  geom_bar(aes(fill = expGroup), stat = "identity", position = "dodge") +
  # geom_text(aes(label = n), angle = 90, hjust = 1.6, position = position_dodge(.9), size = 3) +
  geom_text(aes(x = 1.5, y = 70, label = paste0("p = ", signif(p, 3))), data = pValDF) +
  facet_grid(GRB ~ study + tissue + DE)  + 
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "bottom") + 
  scale_fill_manual(values = COL_RAND) +
  labs(x = "Genes per TAD", y =  "% of TADs")
p
ggsave(p, file = paste0(outPrefix, 
                        ".DEgenes_per_TAD_by_GRB_and_permutation_and_TADsource.LPS_8.pdf"), w = 12, h = 6)



#' 
#' #' Get fisher.test pvalue from a tidy matrix
#' fisherP <- function(df, freq, colx, coly){
#'   
#'   testObj <- df %>% 
#'     select(one_of(c(freq, colx, coly))) %>% 
#'     spread_(key_col = colx, value_col = freq) %>% 
#'     select(-1) %>% 
#'     fisher.test()
#'   
#'   return(testObj$p.value)
#' }
#' 
#' pValDF <- countDF %>%
#'   group_by(GRB, DE) %>% 
#'   do(p = fisherP(., "n", "expGroup", "DE"))

