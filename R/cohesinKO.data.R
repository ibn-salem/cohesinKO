#!/usr/bin/Rscript
#=======================================================================
#
#   parse cohesin KO expression data
#
#=======================================================================

require(gdata)    # to parse .xlsx files (See: http://www.r-bloggers.com/read-excel-files-from-r/)
require(TxDb.Mmusculus.UCSC.mm10.ensGene) # for seqifno object
require(rtracklayer)  # to import bed files as GRanges object with import()
require(tidyverse)  # for tidy data 
require(readxl)   # to read xlsx files
require(stringr)  # to manipulate strings

source("R/GRanges.functions.R")
source("R/screen.grbs.R")

#=======================================================================
# parameters and input data files
#=======================================================================

COHESIN_KO_FILE <- "data/MRC_LMS/RNAseqWTRad21KOMacrophages.RData"
COHESIN_KO_DE_FILE <- "data/MRC_LMS/DifferentialExpr_WTFL_LPS_vs_WTFL_UT_Merged_DESeq2.xlsx"

RUDAN_TAD_FILE <- "data/Rudan2015/mmc2.xlsx"
RAO_MOUSE_TAD_FILE <- "data/Rao2014/GSE63525_CH12-LX_Arrowhead_domainlist.txt.bed.mm9.bed.mm10.bed"
DIXON_MOUSE_TAD_FILES = c(
  Dixon_mESC = "data/Dixon2012/mouse.mESC.mm9.bed.mm10.bed",
  Dixon_cortex = "data/Dixon2012/mouse.cortex.mm9.bed.mm10.bed"
)
GRB_FILE = "data/Harmston2016/mm9_galGal4_70_50.final.bed.mm10.bed"


#=======================================================================
# Parse expression data and reformat
#=======================================================================

# load R file:
cohesinEnv <- new.env()
load(COHESIN_KO_FILE, envir = cohesinEnv)

# convert to tibble and drop non-unique ENSG rows
expDF <- as.tibble(cohesinEnv$FPKM.macrophages) %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE)

# # set ENSG as row names and remove gene name and ID columns
# row.names(expDF) <- expDF[,"ensembl_gene_id"]
# rawExpDF <- expDF[,-(1:2)]

# parse column data and add sample name
cd <- as_tibble(cohesinEnv$colData) %>% 
  mutate(sample = rownames(cohesinEnv$colData)) %>% 
  separate(sample, into = c("rep_name", "condition"), sep="_", remove = FALSE) %>%
  dplyr::select(-rep_name)

expKO <- expDF %>%
  select(1:2, 2 + which(cd$Genotype == "Rad21KO"))

expWT <- expDF %>% 
  select(1:2, 2 + which(cd$Genotype == "WT"))


#=======================================================================
# Parse DE data
#=======================================================================
pairsRow <- read_excel(COHESIN_KO_DE_FILE, n_max = 1, col_names = FALSE)
pairs <- pairsRow %>% 
  t() %>% as_tibble %>% 
  drop_na(V1) %>% unlist()

colRow <- read_excel(COHESIN_KO_DE_FILE, skip = 1, n_max = 1, col_names = FALSE)
cols <- colRow %>% 
  t() %>% as_tibble %>% distinct() %>% 
  filter(! (V1 %in% c("EnsemblID", "GeneSymbol"))) %>% 
  unlist()

colNames <- paste(
  rep(pairs, each = length(cols)),
  rep(cols, length(pairs)),
  sep = "."
)

# parse excel file with differnetial expression data
deDF <- read_excel(COHESIN_KO_DE_FILE, 
                   skip = 2, na = "NA", 
                   col_names = c("EnsemblID", "GeneSymbol", colNames)) 

# data frame with annotation for conditions
conditionDF <- cd %>% 
  select(condition, Genotype, Treatment, Time_hrs) %>% 
  distinct()


# make a tidy data.frame holding for each gene and pairwise condition the
# expression forld change. 
tidyGeneDE <- deDF %>% 
  gather(key = type, value = value, matches(".*_Vs_.*")) %>% 
  separate(type, into = c("comb", "measurement"), sep = "\\.") %>% 
  spread(measurement, value) %>% 
  mutate(DE = factor(
    padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1),
    c(TRUE, FALSE),
    c("DE", "Not DE"))) %>% 
  separate(comb, into = c("cond1", "cond2"), sep = "_Vs_", remove = FALSE) %>% 
  left_join(conditionDF, by = c("cond1" = "condition")) %>%  
  unite(Genotype, Treatment, Time_hrs, col = "cond", sep = ".", remove = FALSE) %>% 
  select(-cond1, -cond2)

dir.create("results", showWarnings = FALSE)
save(deDF, file = "results/deDF.Rdata")
save(tidyGeneDE, file = "results/tidyGeneDE.Rdata")

# load("results/deDE.Rdata")
# load("results/tidyGeneDE.Rdata")

#=======================================================================
# combine replicates
#=======================================================================

# melt data to only one value and one variable column
# add metadata from colData
lexpDF <- expDF %>% 
  gather(3:ncol(expDF), key = "sample", value = "exp") %>% 
  left_join(cd, by = "sample")

repDF <- lexpDF %>% 
  # filter(Genotype == "Rad21KO") %>%
  group_by(Genotype, ensembl_gene_id, Treatment, Time_hrs) %>% 
  summarize(exp_mean = mean(exp, na.rm = TRUE)) %>% 
  unite(Treatment, Time_hrs, col = "condition", sep = "_") %>% 
  spread(condition, exp_mean) %>% 
  left_join(select(expDF, ensembl_gene_id, external_gene_id),
            by = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, external_gene_id, everything())


repKO <- repDF %>% 
  filter(Genotype == "Rad21KO") %>% 
  ungroup() %>% 
  select(-Genotype)


repWT <- repDF %>% 
  filter(Genotype == "WT") %>% 
  ungroup() %>% 
  select(-Genotype)

# # use ENSG as row.names and remove it from columns
# row.names(repKO) <- repKO$ensembl_gene_id
# row.names(repWT) <- repWT$ensembl_gene_id
# 
# repKO$ensembl_gene_id <- NULL
# repWT$ensembl_gene_id <- NULL

expCoDFlist <- list("expWT" = expWT,
                    "expKO" = expKO,
                    "repWT" = repWT,
                    "repKO" = repKO
)

expSource <- tibble(
  name = names(expCoDFlist),
  tissue = rep("Macrophages", 4),
  genotype = c("WT", "Rad21KO", "WT", "Rad21KO"),
  groups = c("treatment.time.replicate", "treatment.time.replicate", "treatment.time", "treatment.time")
)

# test that mean works 
#~ expKO["ENSMUSG00000000028", c("BM1_FL8", "BM2_FL8", "BM3_FL8")]
#~ sum(expKO["ENSMUSG00000000028", c("BM1_FL8", "BM2_FL8", "BM3_FL8")]) / 3
#~ names(repKO)
#~ cohesinEnv$colData
#~ repKO["ENSMUSG00000000028", ]

#===============================================================================
# parse TAD data
#===============================================================================
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
seqInfo <- seqinfo(txdb)

#-------------------------------------------------------------------------------
# Rudan et al. 2015 mouse liver TADs
#-------------------------------------------------------------------------------

# parse bed-like file from .xls table
df <- gdata::read.xls(RUDAN_TAD_FILE, sheet=1)

# convert to GRanges object with 1-based coordinates
Rudan_liver_TAD <-  GRanges(df[,1], IRanges(df[,2]+1, df[,3]), seqinfo=seqInfo)

#-------------------------------------------------------------------------------
# Dixon et 2012  TADs
#-------------------------------------------------------------------------------
DixonTADs <- lapply(DIXON_MOUSE_TAD_FILES, rtracklayer::import.bed, seqinfo=seqInfo)

#-------------------------------------------------------------------------------
# Rudan et al. 2015 mouse liver TADs
#-------------------------------------------------------------------------------
Rao_TADs <- rtracklayer::import.bed(RAO_MOUSE_TAD_FILE, seqinfo = seqInfo)

#-------------------------------------------------------------------------------
# combine all TADs
#-------------------------------------------------------------------------------
allTADs = GRangesList(
  "Rao_CH12-LX"=Rao_TADs,
  "Dixon2012_mESC"=DixonTADs[["Dixon_mESC"]],
  "Dixon2012_cortex"=DixonTADs[["Dixon_cortex"]],
  "VietriRudan_liver"=Rudan_liver_TAD
)

#-------------------------------------------------------------------------------
# GRB-TADs
#-------------------------------------------------------------------------------

# parse GRBs
GRB <- rtracklayer::import.bed(GRB_FILE)

# screen TADs for overlap with GRBs
allTADs <- lapply(allTADs, screen.grbs, GRB)

GRB_TADs <- lapply(allTADs, function(tad) tad[tad$class == "GRB"] )
nonGRB_TADs <- lapply(allTADs, function(tad) tad[tad$class == "nonGRB"] )

# build list with all tads, all GRB TADs and all non-GRB tads
tadList <- c(allTADs, GRB_TADs, nonGRB_TADs)

#-------------------------------------------------------------------------------
# get boundaries of TADs
#-------------------------------------------------------------------------------
boundaryList <- lapply(tadList, getBoundaries)

#-------------------------------------------------------------------------------
# get meta data of TAD sources
#-------------------------------------------------------------------------------

tadSource <- tibble(
  # name = names(allTADs),
  study = c("Rao2014", "Dixon2012", "Dixon2012", "VietriRudan2015"),
  tissue = c("CH12", "mESC", "cortex", "liver")
)

tadSource <- tadSource[rep(1:length(allTADs), 3), ]
tadSource$GRB <- rep(c("all", "GRB", "nonGRB"), each = length(allTADs))

# add unique name to tadList and boundaryList
names(tadList) <- paste0(names(tadList), "_", tadSource$GRB)
names(boundaryList) <- names(tadList)

tadSource <- tadSource %>% 
  mutate(name = names(tadList)) %>% 
  dplyr::select(name, everything())

#-------------------------------------------------------------------------------
# save as .Rdata file
#-------------------------------------------------------------------------------
dir.create("results", showWarnings = FALSE)
save(
  expCoDFlist, 
  expSource, 
  tadList,
  boundaryList,
  tadSource,
  file = "results/Exp_and_TAD_data.Rdata"
  )

