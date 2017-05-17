#!/usr/bin/Rscript
#=======================================================================
#
#   parse cohesin KO expression data
#
#=======================================================================

require(TxDb.Mmusculus.UCSC.mm10.ensGene) # for seqifno object
require(tidyverse)  # for tidy data 
require(DESeq2)     # for differentially expressed genes and fold change

#=======================================================================
# parameters and input data files
#=======================================================================

COHESIN_KO_FILE <- "data/MRC_LMS/RNAseqWTRad21KOMacrophages.RData"



#=======================================================================
# Parse data and reformat
#=======================================================================

# load R file:
cohesinEnv <- new.env()
load(COHESIN_KO_FILE, envir = cohesinEnv)



#=======================================================================
# parse data as SummarizedExperiment
#=======================================================================
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
allGenesGR <- genes(txdb)

mapping <- match(
  cohesinEnv$ReadCounts.macrophages$ensembl_gene_id, 
  allGenesGR$gene_id
)
mapping = mapping[!is.na(mapping)]

# get position in expression data and remove duplicated gene IDs
genesGR <- allGenesGR[mapping]
genesGR <- genesGR[!duplicated(genesGR$gene_id)]

counts <- cohesinEnv$ReadCounts.macrophages %>% 
  filter(ensembl_gene_id %in% names(genesGR)) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
  select(-external_gene_id, -ensembl_gene_id) %>% 
  as.matrix()
colnames(counts) <- NULL

se <- SummarizedExperiment(
  assay = counts,
  colData = cohesinEnv$colData,
  rowRanges = genesGR
)

# relevel the genotype column to have reference (WT) first
se$Genotype <- relevel(se$Genotype, "WT")
# put time as factor
se$Time_hrs <- factor(se$Time_hrs)
se$Replicate <- as.factor(se$Replicate)
se$treatment <- factor(paste(se$Treatment, se$Time_hrs, sep = "_"))
se$group <- as.factor(paste(se$Genotype, se$treatment, sep = "|"))

# # use only LPS treatment
# seLPS <- subset(se, 
#                 select = Treatment %in% c("LPS", "UT") & Time_hrs %in% c(0,8))

#=======================================================================
# Differential expression analysis
#=======================================================================
dds <- DESeqDataSet(se, design = ~ Genotype + treatment)
dds <- DESeq(dds)

resultsNames(dds)

GENOTYPE = "WT"
TREAT = "LPS"
TIME = "8"

fcDF <- tibble(
  gene_id = rowData(se)$gene_id
)

for (GENOTYPE in c("Rad21KO", "WT")) {
  for (TREAT in c("LPS", "LPS+IFNg")) {
    for (TIME in c("2", "8")) {
      
      seSub <- subset(se, 
                      select = Treatment %in% c(TREAT, "UT") & 
                        Time_hrs %in% c(0, TIME) &
                        Genotype == GENOTYPE)
      dds <- DESeqDataSet(seSub, design = ~ Treatment)
      
      dds <- DESeq(dds)
      
      
      res <- results(dds, contrast = c("Treatment", TREAT, "UT"))
      
      # add to gene column    
      resName <- resultsNames(dds)
      colStr <- paste(TREAT, TIME, GENOTYPE, resName, sep = "|")
      
      fcDF <- fcDF %>% 
        mutate(colStr
               
    }
  }
}

res <- results(dds, name = "Time_hrs_8_vs_0")


# rld <- rlog(dds)

