#'
#' Cohesin knock-out analyis. 
#'
#'

require(genepair) # instally by: devtools::install_github("ibn-salem/genepair")
require(biomaRt)  # to download data from ENSMBL
require(TxDb.Mmusculus.UCSC.mm10.ensGene) # for mm10 gene models from ensembl
require(ggplot2)

# load data 
source("R/cohesinKO.data.R")

#===============================================================================
# get gene pairs
#===============================================================================
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
seqInfo <- seqinfo(txdb)

#-------------------------------------------------------------------------------
# Get protein coding genes from ENSEMBL
#-------------------------------------------------------------------------------
allGenesGR <- genes(txdb)

ensemblMouse = useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")

protCodingENSG <- getBM(attributes="ensembl_gene_id", 
                        mart=ensemblMouse, 
                        filters=c("status", "biotype"), 
                        values=list(status="KNOWN", biotype="protein_coding")
                        )

protCodingENSG <- unlist(protCodingENSG)


genesGR <- allGenesGR[allGenesGR$gene_id %in% protCodingENSG]
genesGR <- sort(genesGR)

#-------------------------------------------------------------------------------
# get all close gene pairs
#-------------------------------------------------------------------------------
pairDF <- getAllCisPairs(genesGR, maxDist=10^6, minDist=0)

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
                             filters=c("status", 
                                       "biotype", 
                                       "with_mmusculus_paralog"), 
                             values=list(status="KNOWN", 
                                         biotype="protein_coding",
                                         with_mmusculus_paralog=TRUE),
                             mart=ensemblMouse)

paralogPairs <- data.frame(
  g1 <- match(paralogPairsMouseALL[,1], names(genesGR)),
  g2 <- match(paralogPairsMouseALL[,2], names(genesGR))
)

# which pairs are paralogs
para <- containsGenePairs(pairDF, paralogPairs)

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


# check histogram of distances for all goups
ggplot(pairDF, aes(x=dist, fill=group)) + 
  geom_histogram(breaks=distBreaks) +
  facet_grid(group~., margins=TRUE, scales="free_y") + 
  theme_bw()

#-------------------------------------------------------------------------------
# add TAD annotation
#-------------------------------------------------------------------------------

# get range between tss of gnes as GRanges object
pairGR <- getPairAsGR(pairDF, genesGR)

# iterate over each TAD data set
for(i in 1:length(tadList)){
  
  # get annotation column name
  colname <- paste(as.character(tadSource[i,]), collapse = "|")
  
  # check whether paired genes are in the same TAD
  inTAD <- countOverlaps(pairGR, tadList[[i]], type="within") >= 1
  pairDF[,paste0("TAD_", colname)] <- factor(inTAD, 
                                             c(TRUE, FALSE), 
                                             c("Same TAD", "Not same TAD"))
  
  # check if pairs cross a TAD boundary
  crossBoundary <- countOverlaps(pairGR, boundaryList[[i]]) >= 1
  pairDF[,paste0("Boundary_", colname)] <- factor(
    inTAD, 
    c(TRUE, FALSE), 
    c("Cross Boundary", "Not cross Boundary"))
  
}

#-------------------------------------------------------------------------------
# add expression correlation
#-------------------------------------------------------------------------------

# iterate over each expression data set
for(i in 1:length(expCoDFlist)){
  
  expDF <- expCoDFlist[[i]]
  
  # select matching rows in same order as in genesGR
  subExpDF <- expDF[match(names(genesGR), rownames(expDF)),]
  
  # compute correlation for all pairs
  expCor <- applyToClosePairs(pairDF, genesGR, subExpDF, fun=cor, maxDist=10^6)
  
  dataname <- paste(expSource[i,], collapse="|")
  pairDF[,paste0("expCor_", dataname)] <- expCor
   
}


# create a single (tidy) data frame for all gene pairs with the following columns
# 
# group (paralog, sampled, all)
# dist (abs distance in kb)
# TAD ( sameTAD or notSameTAD)
# Boundary (Cross boundary Not cross boundary)
# tadSource (Study_cell-type)
# expCor (Expression correlation)
# expSource (cell-type_condition_rep..)


