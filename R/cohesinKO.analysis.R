#'
#' Cohesin knock-out downstream analyis. 
#' 
#'

require(ggplot2)  # for plotting
require(RColorBrewer)
require(tidyverse)

# define some parameters
COL_GROUP=brewer.pal(9, "Set1")[c(1,9,2,4)]   # for paralog vs. sampled genes
outPrefix <- "results/v02"

# load data:
load("results/pairDF.Rdata")  # source("cohesinKO.main.R")
load("results/tidyDF.Rdata")  # source("cohesinKO.main.R")


# # separate tadSource and expSource columns
# tidyDF <- tidyDF %>%
#   separate(tadSource, c("study", "tissue", "TADtype"), sep="|", remove=FALSE) %>%
#   transmutate(expSource = gsub("^expCor_", "", expSource)) %>%
#   separate(expSource, c("cell", "genotype", "time", "replicate"), sep="|", remove=FALSE)
  

#-------------------------------------------------------------------------------
# define subset operators for tidyDF
#-------------------------------------------------------------------------------
sExpSource <- tidyDF$expSource == "expCor_Macrophages|WT|treatment.time.replicate"
sTadSource <- tidyDF$tadSource == "Rao2014|CH12|all"

#-------------------------------------------------------------------------------
# count gene pairs by group
#-------------------------------------------------------------------------------

subDF <- subset(tidyDF, sExpSource & sTadSource)

nDF <- subDF %>%
  group_by(group) %>%
  summarise(n=length(group))

p <- ggplot(nDF, aes(x=group, y=n, fill=group)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=n), vjust=-1) +
  theme_bw() + scale_fill_manual(values=COL_GROUP) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position="bottom") +
  xlab("Groups") + ylab("# gene pairs")
ggsave(p, file=paste0(outPrefix, ".n_by_group.hist.pdf"), w=3.5, h=7)

#-------------------------------------------------------------------------------
# compare genomic distance by group
#-------------------------------------------------------------------------------


p <- ggplot(pairDF, aes(x=dist, fill=group)) +
  geom_histogram(breaks=seq(0, 1000, 20), color="black") +
  facet_grid(group~., margins=TRUE, scales="free_y") +
  theme_bw() + scale_fill_manual(values=COL_GROUP) +
  ylab("Gene pairs") + xlab("Genomic distance [kb]")
  # geom_text(aes(x=0, y=0, label=paste0("n=",n)), data=nDF, hjust=-1, vjust=-1) + 
  
ggsave(p, file=paste0(outPrefix, ".dist_by_group.hist.pdf"))


#-------------------------------------------------------------------------------
# pairs in same TAD
#-------------------------------------------------------------------------------
subDF <- subset(tidyDF, sExpSource)
subDF <- as_tibble(subDF)

# separate tadSource and expSource columns
subDF <- subDF %>%
  separate(col=tadSource, into=c("study", "tissue", "TADtype"), remove=FALSE)

# duplicate subDF to add group for all gene pairs 
nSub <- nrow(subDF)
allDF <- rbind(subDF, subDF)

levels(allDF$group) <- c(levels(allDF$group), "all")
allDF$group[seq(nSub+1, 2*nSub)] <- "all"

# pvalDF <- ddply(allDF, .(tad.source), summarize, 
#                 p=fisher.test(table(TAD, group))$p.value)
# 
# pvalDF <- allDF %>%
#   group_by(tadSource) %>%
#   do(fisher.test(table(.TAD, .group))$p.value)


valueDF <- allDF %>%
  group_by(TADtype, study, tissue, group) %>%
  summarise(
    n = sum(TAD == "Same TAD"),
    N = length(TAD),
    percent=100*n/N
    )


p <- ggplot(valueDF, aes(x=group, y=percent)) + 
  geom_bar(aes(fill=group), stat="identity", col="black") + 
  facet_grid(TADtype ~ study+tissue, scales="free_y") +
  scale_fill_manual(values=COL_GROUP, guide_legend(title = "")) +
  theme_bw() + theme(text = element_text(size=15), axis.text.x=element_blank(), legend.position="bottom") + 
  labs(y="Gene pairs in same TAD [%]", x="") + 
  # geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=90), data=pvalDF, size=5) +
  geom_text(aes(label=paste0(round(percent,1),"%" )), vjust=1.5, size=4) + 
  geom_text(aes(label=paste0("n=",n)), vjust=3, size=4)

ggsave(p, file=paste0(outPrefix, ".sameTAD_by_group_and_tadSource.barplot.pdf"), w=14, h=7)


boundaryDF <- allDF %>%
  group_by(TADtype, study, tissue, group) %>%
  summarise(
    n = sum(Boundary == "Not cross Boundary"),
    N = length(Boundary),
    percent=100*n/N
  )


p <- ggplot(boundaryDF, aes(x=group, y=percent)) + 
  geom_bar(aes(fill=group), stat="identity", col="black") + 
  facet_grid(TADtype ~ study+tissue, scales="free_y") +
  scale_fill_manual(values=COL_GROUP, guide_legend(title = "")) +
  theme_bw() + theme(text = element_text(size=15), axis.text.x=element_blank(), legend.position="bottom") + 
  labs(y="Gene pairs not separated by TAD boundary [%]", x="") + 
  # geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=90), data=pvalDF, size=5) +
  geom_text(aes(label=paste0(round(percent,1),"%" )), vjust=1.5, size=4) + 
  geom_text(aes(label=paste0("n=",n)), vjust=3, size=4)
ggsave(p, file=paste0(outPrefix, ".notSeparated_by_group_and_tadSource.barplot.pdf"), w=14, h=7)


#-------------------------------------------------------------------------------
# Expression Correlation
#-------------------------------------------------------------------------------

# separate tadSource and expSource columns
fullDF <- tidyDF %>%
  separate(col=tadSource, into=c("study", "tissue", "TADtype")) %>%
  mutate(expSource = gsub("^expCor_", "", expSource)) %>%
  separate(expSource, c("cell", "genotype", "conditions"), sep="[/|]")

subDF <- subset(fullDF, conditions == "treatment.time.replicate")

#-------------------------------------------------------------------------------
# All gene pairs: cor by Genotype, same TAD , TAD source, and GRB group 
#-------------------------------------------------------------------------------
COL_COMB=brewer.pal(12, "Paired")[c(5,6,7,8,1,2,3,4)]
COL_ALL_GENOTYPE = brewer.pal(12, "Paired")[c(10,9)]

# # calculate p-values
# pvalDF <- ddply(subDF, c("tad.source", "group", tadGroup), summarize, 
#                 p=wilcox.test(expCor ~ expCondition)$p.value)

# #~ nDF <- data.frame(table(subDF[, c("group", "expSource", "TAD")], useNA="ifany"))
# nDF <- ddply(subDF, c("tad.source", "group", tadGroup, "expCondition"), summarize, 
#              n=sum(!is.na(expCor)),
#              nAll=length(expCor),
#              avg=mean(expCor, na.rm=TRUE),
#              sd=sd(expCor, na.rm=TRUE)
# )

p <- ggplot(subDF, aes(x=TAD, y=expCor, color=genotype)) + 
  geom_boxplot(lwd=1.5) + 
  facet_grid(TADtype ~ study+tissue) + 
  theme_bw() + 
  theme(text = element_text(size=15), 
        axis.text.x=element_text(angle = 45, hjust = 1), 
        legend.position="bottom") + 
  scale_color_manual(values=COL_ALL_GENOTYPE, guide_legend(title = "")) +
  guides(fill=guide_legend(title="")) +
  labs(y="Expression Correlation [R]", x="")

  # geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=1.1), data=pvalDF, size=5) +
  # geom_text(aes(label=paste0("n=",n), y=-1), data=nDF) +
  # geom_text(aes(label=paste0("N=",nAll), y=-1.2), data=nDF)

  
ggsave(p, file=paste0(outPrefix, 
                 ".all.expCor_by_genotype_sameTAD_tadSource_GRB.boxplot.pdf"))

#-------------------------------------------------------------------------------
# plot delta of expression values
#-------------------------------------------------------------------------------


# deltaDF <- subDF %>%
#   mutate(gpID = paste0(g1, g2, sep="_")) %>%
#   group_by(gpID, TAD, TADtype, study, tissue) %>%
#   summarise(deltaExp = abs(expCor - expCor) )


