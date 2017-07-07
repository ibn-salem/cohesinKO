#'
#' Cohesin knock-out analyis. 
#'
#'

require(VennDiagram)
require(genepair) # instally by: devtools::install_github("ibn-salem/genepair")
# require(biomaRt)  # to download data from ENSMBL
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

ensemblMouse = biomaRt::useMart(host = "ensembl.org", 
                       biomart = "ENSEMBL_MART_ENSEMBL", 
                       dataset = "mmusculus_gene_ensembl")

protCodingENSG <- biomaRt::getBM(attributes = "ensembl_gene_id", 
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

# save genesGR
save(genesGR, file = "results/genesGR.Rdata")

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
paralogPairsMouseALL = biomaRt::getBM(attributes=paralogParisMouseAttr, 
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
    size = 10000,
    prob = sampWeights,
    replace = FALSE)

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
  inTAD <- countOverlaps(pairGR, tadList[[i]], type = "within") >= 1
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
# generate a tidy data.frame with columns for TAD condition
#-------------------------------------------------------------------------------

# transform pairDF data into a tidy data frame (tidyPairDE) by the follwing opperations
# - combine g1 and g2 to unique gpID column
# - put gpID as first column and remove pralog column
# - combine TAD_*, Boundary_* and nBoundary_* columns
# - separate type (TAD/Boundary/nBoundary) and TAD source
# - spread type (TAD/Boundary/nBoundary) for diffrent sources in single column
# - separate TADsource into study_tissue, study, tissue, and TADtype column
# - remove redundant study_tissue column

tidyPairsTad <- as_tibble(pairDF) %>%
  mutate(gpID = paste(g1, g2, sep = "_")) %>%
  select(gpID, everything(), -paralog) %>% 
  gather(tad_key, tad_value, 
         starts_with("TAD_"), 
         starts_with("Boundary_"), 
         starts_with("nBoundary_")) %>%
  separate(tad_key, 
           into = c("type", "tadSource"), 
           sep = "_", 
           extra = "merge") %>% 
  spread(key = type, value = tad_value) %>%
  separate(
    col = tadSource, 
    into = c("study_tissue", "study", "tissue", "TADtype"), 
    sep = "\\|") %>% 
  mutate(nBoundary = parse_integer(nBoundary)) %>% 
  select(-study_tissue)

# save temporary pairDF data frame
save(tidyPairsTad, file = "results/tidyPairsTad.Rdata")
#load("results/tidyPairsTad.Rdata")

#===============================================================================
# add expression fold changes
#===============================================================================

# select only needed column from tidyGeneDE data.frame
subDE <- tidyGeneDE %>% 
  select(EnsemblID, cond, Genotype, Treatment, Time_hrs, log2FoldChange, DE, reg)

# add ENSG ID to gene paird DF
tidyPairsDE <- as_tibble(pairDF) %>% 
  select(-starts_with("TAD"), -starts_with("Boundary"), -starts_with("nBoundary")) %>% 
  mutate(
    ensg1 = genesGR$gene_id[g1],
    ensg2 = genesGR$gene_id[g2]
  ) %>% 
  select(g1, g2, ensg1, ensg2, everything())

# add fold change for both genes 
tidyPairsDEtmp <- tidyPairsDE %>% 
  left_join(subDE, by = c("ensg1" = "EnsemblID")) %>% 
  left_join(subDE, 
            by = c("ensg2" = "EnsemblID", 
                   "cond" = "cond", 
                   "Genotype" = "Genotype", 
                   "Treatment" = "Treatment", 
                   "Time_hrs" = "Time_hrs"),
            suffix = c("_1", "_2"))

#' Combine two vectors of same factors while ensuring that the samler one is first.
#' 
comb_sorted <- function(x, y, sep = "_"){
  
  # assume factors as input
  stopifnot(is.factor(x) & is.factor(y))
  # assume same set of levels
  stopifnot(all(levels(x) == levels(y)))
  
  firstSmall <- as.ordered(x) <= as.ordered(y)
  
  first <- ifelse(firstSmall, as.character(x), as.character(y))
  second <- ifelse(firstSmall, as.character(y), as.character(x))
  
  stringr::str_c(first, second, sep = sep)
}

# add differnece and log2 ratio to tidyPairsDE for all conditions (combinations)
tidyPairsDE <- tidyPairsDEtmp %>% 
  mutate(
    lfc_diff = abs(log2FoldChange_1 - log2FoldChange_2),
    DE_pair = ifelse(
      is.na(DE_1) | is.na(DE_2),
      NA,
      paste(DE_1, DE_2, sep = "/")),
    DE_pair = ifelse(DE_pair == "Not DE/DE", "DE/Not DE", DE_pair),
    reg_pair = comb_sorted(reg_1, reg_2)
  )

save(tidyPairsDE, file = "results/tidyPairsDE.Rdata")
#load("results/tidyPairsDE.Rdata")

#===============================================================================
# Combine TAD and DE fold change for pairs
#===============================================================================
tidyPairsTadDE <- tidyPairsTad %>%
  left_join(
    select(tidyPairsDE, -(ensg1:group)),
    by = c("g1", "g2")
  )

save(tidyPairsTadDE, file = "results/tidyPairsTadDE.Rdata")
#load("results/tidyPairsTadDE.Rdata")

#-------------------------------------------------------------------------------
# select a single TAD data set and subset to singleTadSourceDF
#-------------------------------------------------------------------------------

# select a single TAD data set
singleTadSource <- tidyPairsTadDE$study == "Rao2014" &
  tidyPairsTadDE$tissue == "CH12" &
  tidyPairsTadDE$TADtype == "all"

singleTadSourceDF <-  subset(tidyPairsTadDE, singleTadSource)
save(singleTadSourceDF, file = "results/singleTadSourceDF.Rdata")

#===============================================================================
# add expression correlation
#===============================================================================

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
  gather(tad_key, tad_value, starts_with("TAD_"), starts_with("Boundary_"), starts_with("nBoundary_")) %>%
  separate(tad_key, into = c("type", "tadSource"), sep = "_", extra = "merge") %>% 
  spread(key = type, value = tad_value) %>% 
  separate(col = tadSource, 
           into = c("study_tissue", "study", "tissue", "TADtype"), 
           sep = "\\|")

save(tmpDF, file = "results/tmpDF.Rdata")

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


#===============================================================================
# get a data frame (TAD2geneDF) that maps each TAD to the containing genes 
#===============================================================================


# get TSS of all genes
tssGR <- resize(genesGR, width = 1, fix = "start", ignore.strand = FALSE)

hits <- map(tadList, findOverlaps, tssGR, ignore.strand = TRUE) %>% 
  map(as.data.frame) %>% 
  map(as_tibble) %>% 
  bind_rows(.id = "name")

# add TAD annotations and gene ID

TAD2geneDF <- hits %>% 
  left_join(tadSource, by = "name") %>% 
  unite(TADid, name, queryHits, sep = "_", remove = FALSE) %>% 
  mutate(ENSG = genesGR$gene_id[subjectHits])

# save TAD2geneDF
save(TAD2geneDF, file = "results/TAD2geneDF.Rdata")
