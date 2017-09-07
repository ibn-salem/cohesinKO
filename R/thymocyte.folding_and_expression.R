#!/usr/bin/Rscript
#===============================================================================
#
#   Analyse thymocyte DKO RNA-seq and Hi-C data
#
#===============================================================================

require(rtracklayer)  # to parse .bed files
require(TxDb.Mmusculus.UCSC.mm9.knownGene) # for seqinfo object
require(biomaRt)      # to get gene coordinates from ensembl
require(tidyverse)    # for gather() function

#===============================================================================
# Some input parametes
#===============================================================================

TAD_FILES <- c(
  "stage_1" = "data/MRC_LMS/thymocyte/TAD_custom/TADs_list_40kb_iced_isHIC_CD69nDP1.bed",
  "stage_4" = "data/MRC_LMS/thymocyte/TAD_custom/TADs_list_40kb_iced_isHIC_CD4SPR1.bed"
)

CD_FILES <- c(
  "cd_stag1_rep1" = "data/MRC_LMS/thymocyte/Hi-C_custom/contact_domain_list_isHIC_CD4SPR1_10kb.txt",
  "cd_stag4_rep1" = "data/MRC_LMS/thymocyte/Hi-C_custom/contact_domain_list_isHIC_CD69nDPR1_10kb.txt",
  "cd_stag4_rep2" = "data/MRC_LMS/thymocyte/Hi-C_custom/contact_domain_list_isHIC_CD69nDPR2_10kb.txt"
)

DE_GENES_FIEL <- "data/MRC_LMS/thymocyte/RNAseq_custom/04_CD69nDPWT1-3_vs_CD69pCD4SPWT1-2.csv"

outPrefix <- "results/thymocyte.TADs_and_DE_by_stage.v02"
dir.create("results", showWarnings = FALSE)

#===============================================================================
# parse data
#===============================================================================

seqInfo <- seqinfo(TxDb.Mmusculus.UCSC.mm9.knownGene)

#-------------------------------------------------------------------------------
# parse TADs
#-------------------------------------------------------------------------------
tadList <- map(TAD_FILES, import.bed) %>% 
  map(function(gr) {
    strand(gr) <- "*"
    return(gr)})

#-------------------------------------------------------------------------------
# parse contact domains
#-------------------------------------------------------------------------------

parseContactDomain <- function(file, ...){
  df <- read_tsv(file, col_types = cols(
    chr1 = col_character(),
    x1 = col_integer(),
    x2 = col_integer(),
    chr2 = col_character(),
    y1 = col_integer(),
    y2 = col_integer(),
    color = col_number(),
    score = col_double(),
    uVarScore = col_double(),
    lVarScore = col_double(),
    upSign = col_double(),
    loSign = col_double()
  ))
  
  GRanges(df$chr1, IRanges(df$x1, df$x2), ...)
}

cdList <- map(CD_FILES, parseContactDomain, seqinfo = seqInfo)

# combine TADs and contact dominas
domainList <- c(tadList, cdList)

domain_meta <- tibble(
  name = c(names(tadList), names(cdList)),
  stage = c(1, 4, 1, 4, 4),
  type = c("TAD", "TAD", rep("contact_domain", 3)),
  HiC_rep = c(1, 1, 1, 1, 2)
)

#-------------------------------------------------------------------------------
# parse expression data
#-------------------------------------------------------------------------------
expDat <- read_csv(DE_GENES_FIEL, col_types = cols(
  ensembl_gene_id = col_character(),
  baseMean = col_double(),
  log2FoldChange = col_double(),
  lfcSE = col_double(),
  stat = col_double(),
  pvalue = col_double(),
  padj = col_double(),
  mgi_symbol = col_character()
)) 

# remove non-unique gene ids
expDat <- expDat %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE)

# add label for regulated genes
expDat <- expDat %>% 
  mutate(
    DE = factor(
      padj <= 0.05 & (log2FoldChange > 1 | log2FoldChange < -1),
      c(TRUE, FALSE),
      c("DE", "Not DE")),
    reg = factor(
      ifelse(DE == "DE", 
             ifelse(log2FoldChange > 1, "up", "down"),
             "unchanged"),
      c("up", "down", "unchanged"))
  )
    


#-------------------------------------------------------------------------------
# get genes from ENSEMBL for mm9
#-------------------------------------------------------------------------------
# Ensmbl archive with mm9 genome: http://may2012.archive.ensembl.org
ensemblmm9 <- useMart(host = "may2012.archive.ensembl.org",
                         biomart = "ENSEMBL_MART_ENSEMBL", 
                         dataset = "mmusculus_gene_ensembl", 
                         verbose = FALSE)

geneAttributes = c("ensembl_gene_id", "external_gene_id", "chromosome_name", 
                   "start_position", "end_position", "strand", "status", 
                   "gene_biotype")

genes = as_tibble(getBM(
  attributes = geneAttributes, 
  mart = ensemblmm9, 
  filters = list("ensembl_gene_id" = expDat$ensembl_gene_id)
))


# make GRanges object for tss of all known prot coding genes
tssCoord <- ifelse(genes$strand == 1, genes$start_position, genes$end_position)

tssGR <- GRanges(paste0("chr", genes$chromosome_name),
                 IRanges(tssCoord, tssCoord),
                 strand = ifelse(genes$strand == 1, '+', '-'), 
                 dplyr::select(genes, -chromosome_name, -start_position, -end_position, -strand),
                 seqinfo = seqInfo
)

#===============================================================================
# Analyse gene expression in TADs
#===============================================================================

tssDF <- as_tibble(as.data.frame(mcols(tssGR)))
# add column indicating if tss of gene is within TAD


olDF <- as_tibble(
  map(domainList, function(gr) countOverlaps(tssGR, gr) > 0)
)

tssDF <- tssDF %>% 
  bind_cols(olDF)


# add expression data

realExp <- expDat %>% 
  select( ensembl_gene_id, DE) %>% 
  mutate(
    group = "DiffExp",
    rand_rep = 1
    )


# iterate for randomizations
randDF <- bind_rows(map(1:10, function(i){

  randExp <- expDat %>% 
    select( ensembl_gene_id, DE) %>% 
    mutate(
      DE = sample(DE),
      group = "random",
      rand_rep = i
    )
}))

tssDF <- tssDF %>% 
  left_join(bind_rows(realExp, randDF), by = "ensembl_gene_id")


# combine differnt TAD/CD data sets
tssDeDomainDF <- tssDF %>% 
  gather(names(domainList), key = "TAD_type", value = "inTAD")

save(tssDeDomainDF, file = paste0(outPrefix, "tssDeDomainDF.Rdata"))

#-------------------------------------------------------------------------------
# count number of genes in TADs
#-------------------------------------------------------------------------------

countDF <- tssDeDomainDF %>% 
  count(group, rand_rep, TAD_type, DE, inTAD) %>% 
  mutate(
    percent = 100 * n / n_distinct(tssDeDomainDF$ensembl_gene_id)
  ) 

# combine randomization replicates by mean and sd  
combinedDF <- countDF %>% 
  group_by(group, TAD_type, DE, inTAD) %>% 
  summarize(
    n_mean = mean(n, na.rm = TRUE),
    n_sd = sd(n, na.rm = TRUE),
    percent_mean = mean(percent, na.rm = TRUE),
    percent_sd = sd(percent, na.rm = TRUE)
  ) %>% 
  left_join(domain_meta, by = c("TAD_type" = "name")) %>% 
  mutate(stage = as.factor(stage))

write_tsv(combinedDF, paste0(outPrefix, ".DEgeens_inDomain.values.tsv"))

inDomainDF <- combinedDF %>% 
  filter(inTAD) %>% 
  filter(DE == "DE")

# percent
p <- ggplot(inDomainDF, aes(x = stage, y = percent_mean, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = signif(percent_mean, 2)), 
            position = position_dodge(width = .9),
            vjust = -.6) + 
  geom_errorbar(aes(ymin = percent_mean - percent_sd, ymax = percent_mean + percent_sd),
                size = .3,    # Thinner lines
                width = .2,
                position = position_dodge(.9)) +
  facet_grid(. ~ type + HiC_rep, scales = "free_x") + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set1") + 
  labs(y = "Percent DE genes in domain")

ggsave(p, filename = paste0(outPrefix, ".DEgeens_inDomain.percent.barplot.pdf"), w = 5, h = 5)

# gene counts
p <- ggplot(inDomainDF, aes(x = stage, y = n_mean, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = n_mean), 
            position = position_dodge(width = .9),
            vjust = -.6) + 
  geom_errorbar(aes(ymin = n_mean - n_sd, ymax = n_mean + n_sd),
                size = .3,    # Thinner lines
                width = .2,
                position = position_dodge(.9)) +
  facet_grid(. ~ type + HiC_rep, scales = "free_x") + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set1") + 
  labs(y = "DE genes in domain")

ggsave(p, filename = paste0(outPrefix, ".DEgeens_inDomain.n.barplot.pdf"), w = 5, h = 5)

#===============================================================================
# Analyse overalp of TADs
#===============================================================================

# get numer of TADs without any overlap
unique_TAD_1 <- sum(countOverlaps(tadList[[1]], tadList[[2]], ignore.strand = TRUE) == 0)
unique_TAD_4 <- sum(countOverlaps(tadList[[2]], tadList[[1]], ignore.strand = TRUE) == 0)

# compute reciprocal overlap of TADs
recipTADhits <- reciprocalOverlaps(tadList[[1]], tadList[[2]], 0.5, ignore.strand = TRUE)
unique50_TAD_1 = sum(countLnodeHits(recipTADhits) == 0)
unique50_TAD_4 = sum(countRnodeHits(recipTADhits) == 0)


# same for contact domains
cd1 <- domainList[["cd_stag1_rep1"]]
cd2 <- domainList[["cd_stag4_rep1"]]

unique_CD_1 <- sum(countOverlaps(cd1, cd2, ignore.strand = TRUE) == 0)
unique_CD_4 <- sum(countOverlaps(cd2, cd1, ignore.strand = TRUE) == 0)

# compute reciprocal overlap of TADs
recipCDhits <- reciprocalOverlaps(cd1, cd2, 0.5, ignore.strand = TRUE)
unique50_CD_1 = sum(countLnodeHits(recipCDhits) == 0)
unique50_CD_4 = sum(countRnodeHits(recipCDhits) == 0)


# ol <- findOverlaps(tadList[[1]], tadList[[2]], ignore.strand = TRUE)

tad_stats_DF <- tibble(
  name = names(tadList),
  n =  map_int(tadList, length),
  n_uniq_TAD = c(unique_TAD_1, unique_TAD_4),
  n_uniq50_TAD = c(unique50_TAD_1, unique50_TAD_4),
  width_median = map_dbl(tadList, function(gr) median(width(gr))),
  width_mean = map_dbl(tadList, function(gr) mean(width(gr))),
  sd_median = map_dbl(tadList, function(gr) sd(width(gr))),
  genome_cov = map_dbl(tadList, function(gr) 
    sum(as.numeric(width(GenomicRanges::intersect(gr, GRanges(seqInfo), ignore.strand = TRUE))))
  ),
  genome_cov_percent = 100 * genome_cov / sum(as.numeric(width(GRanges(seqInfo))))
)

write_tsv(tad_stats_DF, paste(outPrefix, ".tad_stats.tsv"))

#-------------------------------------------------------------------------------
# Group genes by beeing in stag1 TAD but not in stage2
#-------------------------------------------------------------------------------

# for each gene check if it is inside TAD in both, no or only one stage
tssTADchangeDF <- tssDF %>% 
  mutate(
    TAD = case_when(
      stage_1 & ! stage_4 ~ "stage_1_only",
      stage_4 & ! stage_1 ~ "stage_4_only",
      stage_1 & stage_4 ~ "both_inside",
      !stage_1 & !stage_4 ~ "both_outside"
    ),
    CD = case_when(
      cd_stag1_rep1 & ! cd_stag4_rep1 ~ "stage_1_only",
      cd_stag4_rep1 & ! cd_stag1_rep1 ~ "stage_4_only",
      cd_stag1_rep1 & cd_stag4_rep1 ~ "both_inside",
      !cd_stag1_rep1 & ! cd_stag4_rep1 ~ "both_outside"
    )) %>% 
  mutate(
    TAD = factor(TAD, c("stage_1_only", "stage_4_only", "both_inside", "both_outside")),
    CD = factor(CD, c("stage_1_only", "stage_4_only", "both_inside", "both_outside"))
  ) %>% 
  gather(key = "TAD_type", value = "change", TAD, CD)

# count number of DE genes and percent for each group
changeCountDF <- tssTADchangeDF %>% 
  group_by(group, rand_rep, TAD_type, change, DE) %>% 
  summarize(
    n = n()
  ) %>% 
  mutate(
    percent = n / sum(n) * 100
  ) 

# combine counts by replicates
combinedChangeDF <- changeCountDF %>% 
  group_by(group, TAD_type, change, DE) %>% 
  summarize(
    n_mean = mean(n, na.rm = TRUE),
    n_sd = sd(n, na.rm = TRUE),
    percent_mean = mean(percent, na.rm = TRUE),
    percent_sd = sd(percent, na.rm = TRUE)
  ) 

write_tsv(combinedChangeDF, paste(outPrefix, ".DEgeens_changed_TADs.tsv"))

# percent
p <- ggplot(combinedChangeDF, aes(x = change, y = percent_mean, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = signif(percent_mean, 2)), 
            position = position_dodge(width = .9),
            vjust = -.6) + 
  geom_errorbar(aes(ymin = percent_mean - percent_sd, ymax = percent_mean + percent_sd),
                size = .3,    # Thinner lines
                width = .2,
                position = position_dodge(.9)) +
  facet_grid(TAD_type ~ DE, scales = "free_y") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1") + 
  labs(y = "Percent of genes")
ggsave(p, filename = paste0(outPrefix, ".DEgeens_changed_TADs.percent.barplot.pdf"), w = 6, h = 6)

# percent
p <- ggplot(combinedChangeDF, aes(x = change, y = n_mean, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = signif(n_mean, 2)), 
            position = position_dodge(width = .9),
            angle = 90, hjust = "inward") + 
  geom_errorbar(aes(ymin = n_mean - n_sd, ymax = n_mean + n_sd),
                size = .3,    # Thinner lines
                width = .2,
                position = position_dodge(.9)) +
  facet_grid(TAD_type ~ DE, scales = "free_y") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1") + 
  labs(y = "Number of genes")
ggsave(p, filename = paste0(outPrefix, ".DEgeens_changed_TADs.n.barplot.pdf"), w = 6, h = 6)

