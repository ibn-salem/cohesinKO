#'
#' Cohesin knock-out analyis. 
#'
#'

require(genepair) # instally by: devtools::install_github("ibn-salem/genepair")
require(biomaRt)  # to download data from ENSMBL
require(TxDb.Mmusculus.UCSC.mm10.ensGene) # for mm10 gene models from ensembl
require(tidyverse)    # for gather() function

# load data 
# source("R/cohesinKO.data.R")
load("results/Exp_and_TAD_data.Rdata")
load("results/tidyGeneDE.Rdata")

#===============================================================================
# get gene pairs
#===============================================================================
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
seqInfo <- seqinfo(txdb)

#-------------------------------------------------------------------------------
# Get protein coding genes from ENSEMBL
#-------------------------------------------------------------------------------
allGenesGR <- genes(txdb)

expGenes <- unique(tidyGeneDE$EnsemblID)

ensemblMouse = useMart(host = "ensembl.org", 
                       biomart = "ENSEMBL_MART_ENSEMBL", 
                       dataset = "mmusculus_gene_ensembl")

protCodingENSG <- getBM(attributes = "ensembl_gene_id", 
                        mart = ensemblMouse, 
                        filters = c("status", "biotype"), 
                        values = list(status = "KNOWN", 
                                      biotype = "protein_coding")
)

protCodingENSG <- unlist(protCodingENSG)

#-------------------------------------------------------------------------------
# analyse overlap of Gene IDs in different data sets----------------------------
#-------------------------------------------------------------------------------
vennList <- list(
  expression = unique(tidyGeneDE$EnsemblID),
  TxDb.Mmusculus.UCSC.mm10.ensGene = names(allGenesGR),
  protCodingENSG = protCodingENSG
)

require(VennDiagram)
venn.diagram(
  x = vennList, 
  filename = "results/ENSG_overlap.venn.png",
  imagetype = "png",
  euler.d = TRUE, scaled = TRUE,
  lwd = 3,
  lty = 'blank',
  fill = c("cornflowerblue", "green", "yellow"),
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 2,
  cat.fontfamily = "sans",
  rotation = 1
)


#-------------------------------------------------------------------------------
# get only those genes that are in the expression data set and annotated as protein coding in ensembl
#-------------------------------------------------------------------------------
genesGR <- allGenesGR[
  allGenesGR$gene_id %in% protCodingENSG & 
    allGenesGR$gene_id %in% expGenes]

genesGR <- sort(genesGR)

#-------------------------------------------------------------------------------
# get all close gene pairs
#-------------------------------------------------------------------------------
pairDF <- getAllCisPairs(genesGR, maxDist = 10^6, minDist = 0)

# convert distance in absolut kb distances
pairDF$dist <- abs(pairDF$dist) / 1000

#-------------------------------------------------------------------------------
# Mouse paralog gene pairs
#-------------------------------------------------------------------------------
paralogParisMouseAttr = c("ensembl_gene_id", 
                          "mmusculus_paralog_ensembl_gene", 
                          "mmusculus_paralog_perc_id", 
                          "mmusculus_paralog_perc_id_r1", 
                          "mmusculus_paralog_dn", 
                          "mmusculus_paralog_ds")     

# download all paralog pairs
paralogPairsMouseALL = getBM(attributes=paralogParisMouseAttr, 
                             filters = c("status", 
                                       "biotype", 
                                       "with_mmusculus_paralog"), 
                             values = list(status = "KNOWN", 
                                         biotype = "protein_coding",
                                         with_mmusculus_paralog = TRUE),
                             mart = ensemblMouse)

paralogPairs <- data.frame(
  g1 = match(paralogPairsMouseALL[,1], names(genesGR)),
  g2 = match(paralogPairsMouseALL[,2], names(genesGR))
)

# which pairs are paralogs
pairDF[,"paralog"] <- containsGenePairs(pairDF, paralogPairs)

#-------------------------------------------------------------------------------
# add sampled pairs with same distance distribution as paralogs
#-------------------------------------------------------------------------------

# get distance for paralog and non-paralog gene pairs
distPara <- pairDF[pairDF$paralog, "dist"]
distNonPara <- pairDF[!pairDF$paralog, "dist"]

# devide distance range in 50 bins a 20kb
distBreaks <- seq(0, 1000, 20)

# get samplig weights according to observed distance in paralogs
sampWeights <- weightsByBin(distPara, distNonPara, distBreaks)

# sample from non-paralog pairs with weights
sampIdx <- sample(
    which(!pairDF$paralog), 
    size=10000,
    prob=sampWeights,
    replace=FALSE)

# add group column to pairDF data.frame
pairDF[,"group"] <- ifelse(pairDF[,"paralog"], "paralog", "non-paralog")
pairDF[sampIdx, "group"] <- "sampled"
pairDF[,"group"] <- factor(pairDF[,"group"], c("paralog", "sampled", "non-paralog"))


# # check histogram of distances for all goups
# ggplot(pairDF, aes(x=dist, fill=group)) + 
#   geom_histogram(breaks=distBreaks) +
#   facet_grid(group~., margins=TRUE, scales="free_y") + 
#   theme_bw()

#-------------------------------------------------------------------------------
# add TAD annotation
#-------------------------------------------------------------------------------

# get range between tss of gnes as GRanges object
pairGR <- getPairAsGR(pairDF, genesGR)

# iterate over each TAD data set
for (i in 1:length(tadList)) {
  
  # get annotation column name
  colname <- paste(as.character(tadSource[i,]), collapse = "|")
  
  # check whether paired genes are in the same TAD
  inTAD <- countOverlaps(pairGR, tadList[[i]], type="within") >= 1
  pairDF[,paste0("TAD_", colname)] <- factor(inTAD, 
                                             c(TRUE, FALSE), 
                                             c("Same TAD", "Not same TAD"))
  
  # check if pairs cross a TAD boundary
  crossBoundary <- countOverlaps(pairGR, boundaryList[[i]])
  pairDF[,paste0("Boundary_", colname)] <- factor(
    crossBoundary > 1, 
    c(TRUE, FALSE), 
    c("Cross Boundary", "Not cross Boundary"))
  
  # number of boundaries between pairs
  pairDF[,paste0("nBoundary_", colname)] <- crossBoundary
}

# save temporary pairDF data frame
save(pairDF, genesGR, file = "results/pairDF_genesGR_TMP_TAD.Rdata")
#load("results/pairDF_genesGR_TMP_TAD.Rdata")

#-------------------------------------------------------------------------------
# add expression fold changes
#-------------------------------------------------------------------------------

# select only needed column from tidyGeneDE data.frame
subDE <- tidyGeneDE %>% 
  select(EnsemblID, cond, log2FoldChange)

# add ENSG ID to gene paird DF
pairDF <- as_tibble(pairDF) %>% 
  mutate(
    ensg1 = genesGR$gene_id[g1],
    ensg2 = genesGR$gene_id[g2]
  ) %>% 
  select(g1, g2, ensg1, ensg2, everything())

# add fold change for both genes 
pairDF <- pairDF %>% 
  left_join(subDE, by = c("ensg1" = "EnsemblID")) %>% 
  left_join(subDE, 
            by = c("ensg2" = "EnsemblID", "cond" = "cond"),
            suffix = c("_1", "_2"))

# add differnece and log2 ratio to pairDF for all conditions (combinations)
pairDF <- pairDF %>% 
  mutate(lfc_diff = abs(log2FoldChange_1 - log2FoldChange_2)) %>% 
  mutate(lfc_lfc = log2(log2FoldChange_1 / log2FoldChange_2))

save(pairDF, genesGR, file = "results/pairDF_genesGR_TMP_DE.Rdata")
#load("results/pairDF_genesGR_TMP_DE.Rdata")

#--------------------

subPairDF <- pairDF %>% 
  filter(cond == "Rad21KO.LPS.2")

subPairDF %>% ggplot(aes(x = lfc_diff)) + 
  geom_histogram()

subPairDF %>% ggplot(aes(x = lfc_lfc)) + 
  geom_histogram()

subPairDF %>% ggplot(aes(x = log2FoldChange_1)) + 
  geom_histogram()

subPairDF %>% 
  ggplot(aes(x = log2FoldChange_1, y = log2FoldChange_2)) + 
  geom_point(alpha = 0.2)



subPairDF <- pairDF %>% 
  filter(cond %in% c("Rad21KO.LPS.2", "WT.LPS.2"))

subPairDF %>% 
  ggplot(aes(color = cond, x = cond, y = lfc_diff)) + 
  facet_grid(. ~ `Boundary_Rao2014|CH12|GRB`) + 
  geom_boxplot()
  # stat_ecdf(geom = "step", pad = FALSE)

subPairDF %>% 
  ggplot(aes(color = cond, x = cond, y = lfc_lfc)) + 
  facet_grid(. ~ `Boundary_Rao2014|CH12|GRB`) + 
  geom_boxplot()
  # stat_ecdf(geom = "step", pad = FALSE)



#--------------------



#-------------------------------------------------------------------------------
# generate a tidy data.frame with columns for TAD condition and DE conditions
#-------------------------------------------------------------------------------

# transform pairDF data into a tidy data frame (tidyPairDE) by the follwing opperations
# - combine g1 and g2 to unique gpID column
# - put gpID as first column and remove pralog and single gID columns
# - combine TAD_ and Boundary_ column
# - separate type (TAD/Boundary) and TAD source
# - spread type (TAD/Boundary) for diffrent sources in single column
# - separate TADsource into study, tissue, and TADtype column
# - combine expCor_ and expDiff columns
# - separated (expCor/expDiff) into type and tadSource
# - spread by (expCor/expDiff)
# - separate expSource into cell, genotype, and condtion column
# - remove redundant "expCor" prefix column

tmpDF <- pairDF %>%
  mutate(gpID = paste(g1, g2, sep = "_")) %>%
  select(gpID, everything(), -paralog, -g1, -g2) %>% 
  gather(tad_key, tad_value, starts_with("TAD_"), starts_with("Boundary_")) %>%
  separate(tad_key, into = c("type", "tadSource"), sep = "_", extra = "merge") %>% 
  spread(key = type, value = tad_value) %>%
  separate(col = tadSource, into = c("study", "tissue", "TADtype"))


#-------------------------------------------------------------------------------
# add expression correlation
#-------------------------------------------------------------------------------

# get ordered list of gene ids
genesDF <- tibble(id = genesGR$gene_id)

# convert pairDF into a tibble for easiers handling and printig
pairDF <- as_tibble(pairDF)

# use only the replicate averaged expression data
expCoDFlist = expCoDFlist[c("repWT", "repKO")]
expSource = expSource %>% filter(name %in% c("repWT", "repKO"))

# iterate over each expression data set
for (i in 1:length(expCoDFlist)) {
  
  expDF <- expCoDFlist[[i]]
  
  # select matching rows in same order as in genesGR
  subExpDF <- expDF %>% 
    right_join(genesDF, by = c("ensembl_gene_id" = "id")) %>% 
    select(-external_gene_id, -ensembl_gene_id) %>% 
    as.data.frame()
  
  # compute correlation for all pairs
  expCor <- applyToClosePairs(pairDF, genesGR, subExpDF, fun = cor, maxDist = 10^6)
  
  dataname <- paste(expSource[i, ], collapse = "|")
  pairDF[,paste0("expCor_", dataname)] <- expCor
   
}

#' coefficient of variation
#' 
#' Computes coefficient of variation (CV), also known as relative standard
#' deviation (RSD). See
#' \url{https://en.wikipedia.org/wiki/Coefficient_of_variation}.
#' 
#' @param x numeric vector
cv <- function(x, na.rm = FALSE){
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm) * 100
}


# iterate over each expression data set
for (i in 1:length(expCoDFlist)){
  
  expDF <- expCoDFlist[[i]]

  # select matching rows in same order as in genesGR
  subExpDF <- expDF %>% 
    right_join(genesDF, by = c("ensembl_gene_id" = "id"))
  
  tidyExpDF <- subExpDF %>% 
    gather(key = sample, value = exp, -external_gene_id, -ensembl_gene_id)
  
  genesSummaryDF <- tidyExpDF %>% 
    group_by(ensembl_gene_id) %>% 
    summarize(
      n = n(),
      mean_exp = mean(exp, na.rm = TRUE),
      sd_exp = sd(exp, na.rm = TRUE),
      cv = cv(exp, na.rm = TRUE)
    ) %>% 
    mutate(id = match(ensembl_gene_id, genesDF$id))
  
  # add cv of gene 1 and 2 in pair and abs diff
  newPairDF <- select(pairDF, g1, g2) %>% 
    left_join(select(genesSummaryDF,
                     id,
                     g1_mean_exp = mean_exp,
                     g1_sd_exp = sd_exp,
                     g1_cv = cv), by = c("g1" = "id")) %>% 
    left_join(select(genesSummaryDF,
                     id,
                     g2_mean_exp = mean_exp,
                     g2_sd_exp = sd_exp,
                     g2_cv = cv), by = c("g2" = "id")) %>% 
    mutate(cvDiff = abs(g2_cv - g1_cv))
  
  dataname <- paste(expSource[i,], collapse = "|")
  pairDF[,paste0("cvDiff_", dataname)] <- newPairDF$cvDiff
  
}

#-------------------------------------------------------------------------------
# save pairDF 
#-------------------------------------------------------------------------------
save(pairDF, file = "results/pairDF.Rdata")
#load("results/pairDF.Rdata")

#-------------------------------------------------------------------------------
# create a single (tidy) data frame for all pairs with the following columns
#-------------------------------------------------------------------------------

# Aim for the folowing columns:
#
# group (paralog, sampled, all)
# dist (abs distance in kb)
# TAD ( sameTAD or notSameTAD)
# Boundary (Cross boundary Not cross boundary)
# tadSource (Study_cell-type)
# expCor (Expression correlation)
# expSource (cell-type_condition_rep..)

# use tidyr::gather() function as nicely descibed here:
# http://r4ds.had.co.nz/tidy-data.html#gathering


# transform pairDF data into tidyDF by the follwing opperations
# - combine g1 and g2 to unique gpID column
# - put gpID as first column and remove pralog and single gID columns
# - combine TAD_ and Boundary_ column
# - separate type (TAD/Boundary) and TAD source
# - spread type (TAD/Boundary) for diffrent sources in single column
# - separate TADsource into study, tissue, and TADtype column
# - combine expCor_ and expDiff columns
# - separated (expCor/expDiff) into type and tadSource
# - spread by (expCor/expDiff)
# - separate expSource into cell, genotype, and condtion column
# - remove redundant "expCor" prefix column

tmpDF <- pairDF %>%
  mutate(gpID = paste(g1, g2, sep = "_")) %>%
  select(gpID, everything(), -paralog, -g1, -g2) %>% 
  gather(tad_key, tad_value, starts_with("TAD_"), starts_with("Boundary_")) %>%
  separate(tad_key, into = c("type", "tadSource"), sep = "_", extra = "merge") %>% 
  spread(key = type, value = tad_value) %>%
  separate(col = tadSource, into = c("study", "tissue", "TADtype"))
  
tidyDF <- tmpDF  %>% 
  gather(exp_key, exp_value, matches("expCor_.|cvDiff_.")) %>%
  separate(col = exp_key, into = c("exp_type", "expSource"), sep = "_") %>% 
  spread(key = exp_type, value = exp_value) %>%
  separate(expSource, into = c("DFname", "cell", "genotype", "conditions"), sep = "[|]") %>% 
  select(-DFname)

  
  # gather(key = expSource, value = expCor, matches("expCor_.")) %>%
  # mutate(expSource = gsub("^expCor_", "", expSource)) %>%
  # separate(expSource, c("cell", "genotype", "conditions"), sep = "[/|]")

# save tidyDF
save(tidyDF, file = "results/tidyDF.Rdata")


