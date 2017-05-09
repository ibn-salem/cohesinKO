#'
#' Cohesin knock-out downstream analyis. 
#' 
#'

require(RColorBrewer) # for picking nice colors
require(tidyverse)    # to format, sumarize and plot tidy data

source("R/colors.functions.R")

COL_GROUP = brewer.pal(9, "Set1")[c(1,9,2,4)]   # for paralog vs. sampled genes
COL_COMB = brewer.pal(12, "Paired")[c(5,6,7,8,1,2,3,4)]
pairedCols <- brewer.pal(12, "Paired")
COL_GROUP_PAIRED <- c(pairedCols[c(6,5)], "azure4", "azure3", pairedCols[c(2,1,  10,9)])
COL_ALL_GENOTYPE = brewer.pal(12, "Paired")[c(10,9)]
COL_TAD_GENOTYPE = brewer.pal(12, "Paired")[c(8,7,  4,3)] # orange, green

# CONDITIONS="treatment.time.replicate"
CONDITIONS="treatment.time"

# define some parameters
outPrefix <- ifelse(CONDITIONS == "treatment.time", 
                    "results/v03.avgRep", 
                    "results/v03")


# load data:
#load("results/pairDF.Rdata")  # source("cohesinKO.main.R")
load("results/tidyDF.Rdata")  # source("cohesinKO.main.R")

#===============================================================================
# define subset operators for tidyDF
#===============================================================================

# select a single expression data set
singleExpSource <- tidyDF$cell == "Macrophages" &
  tidyDF$genotype == "WT" &
  tidyDF$conditions == CONDITIONS

# select a single TAD data set
singleTadSource <- tidyDF$study == "Dixon2012" &
  tidyDF$tissue == "mESC" &
  tidyDF$TADtype == "all"

# check that each pair is included once
# sum(singleExpSource & singleTadSource) == nrow(pairDF)

#===============================================================================
# count gene pairs by group
#===============================================================================

subDF <- subset(tidyDF, singleExpSource & singleTadSource)

nDF <- subDF %>%
  group_by(group) %>%
  summarise(n=length(group))

p <- ggplot(nDF, aes(x=group, y=n, fill=group)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=n), vjust=-1) +
  theme_bw() + scale_fill_manual(values=COL_GROUP) +
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position="none") +
  xlab("Groups") + ylab("# gene pairs")
ggsave(p, file=paste0(outPrefix, ".n_by_group.hist.pdf"), w = 3.5, h = 7)

#===============================================================================
# compare genomic distance by group
#===============================================================================

p <- ggplot(subDF, aes(x=dist, fill=group)) +
  geom_histogram(breaks=seq(0, 1000, 20), color="black") +
  facet_grid(group~., margins=TRUE, scales="free_y") +
  theme_bw() + scale_fill_manual(values=COL_GROUP, guide=FALSE) +
  theme(text = element_text(size=15)) + 
  ylab("Gene pairs") + xlab("Genomic distance [kb]")
  
ggsave(p, file=paste0(outPrefix, ".dist_by_group.hist.pdf"))

#===============================================================================
# pairs in same TAD
#===============================================================================
# use multiple TAD data sets, but only a single expression data set
subDF <- subset(tidyDF, singleExpSource)
#subDF <- as_tibble(subDF)

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

sameTadDF <- allDF %>%
  group_by(TADtype, study, tissue, group) %>%
  summarise(
    n = sum(TAD == "Same TAD"),
    N = length(TAD),
    percent=100*n/N
    )


p <- ggplot(sameTadDF, aes(x=group, y=percent)) + 
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
    percent = 100*n/N
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


#===============================================================================
# Expression Correlation
#===============================================================================

#-------------------------------------------------------------------------------
# All gene pairs groups: cor by Genotype and same TAD 
#-------------------------------------------------------------------------------

# select only correlation across replicates
subDF <- subset(tidyDF, singleTadSource & 
                  tidyDF$conditions == CONDITIONS)

# calculate p-values for differnece in genotype
pvalDF <- subDF %>%
  group_by(TAD) %>%
  do(w = wilcox.test(expCor ~ genotype, data=.)) %>%
  summarise(TAD, p=w$p.value)
  
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


ggsave(p, file=paste0(outPrefix, 
                      ".all_groups.expCor_by_genotype_sameTAD_Dixon_mESC.boxplot.pdf"), w=3.5, h=7)

#-------------------------------------------------------------------------------
# All gene pairs groups: cor by Genotype and not cross Boundary 
#-------------------------------------------------------------------------------

# calculate p-values for differnece in genotype
pvalDF <- subDF %>%
  group_by(Boundary) %>%
  do(w = wilcox.test(expCor ~ genotype, data=.)) %>%
  summarise(Boundary, p = w$p.value)

p <- ggplot(subDF, aes(x=Boundary, y=expCor, color=genotype)) + 
  geom_boxplot(lwd=1.5) + 
  facet_grid(TADtype ~ study+tissue) + 
  theme_bw() + 
  theme(text = element_text(size=15), 
        axis.text.x=element_text(angle = 45, hjust = 1), 
        legend.position="bottom") + 
  scale_color_manual(values=COL_ALL_GENOTYPE, guide_legend(title = "")) +
  guides(fill=guide_legend(title="")) +
  labs(y="Expression Correlation [R]", x="")

ggsave(p, file=paste0(outPrefix, 
                      ".all_groups.expCor_by_genotype_notSeparated_Dixon_mESC.boxplot.pdf"), w=3.5, h=7)


#-------------------------------------------------------------------------------
# Cor by Genotype, same TAD, and genepair groups 
#-------------------------------------------------------------------------------

# select only correlation across replicates
subDF <- subset(tidyDF, singleTadSource & 
                  tidyDF$conditions == CONDITIONS)

# duplicate to add "all" group
nSub <- nrow(subDF)
allDF <- rbind(subDF, subDF)
levels(allDF$group) <- c(levels(allDF$group), "all")
allDF$group[seq(nSub+1, 2*nSub)] <- "all"

allDF <- allDF %>%
  arrange(group, genotype) %>%
  mutate(genotype_group = paste(group, genotype, sep="_"))  %>%
  mutate(genotype_group = factor(genotype_group, unique(genotype_group)))
  

p <- ggplot(allDF, aes(x=TAD, y=expCor, color=genotype_group)) + 
  geom_boxplot(lwd=1.5) + 
  facet_grid(TADtype ~ group) + 
  theme_bw() + 
  theme(text = element_text(size=15), 
        axis.text.x=element_text(angle = 45, hjust = 1), 
        legend.position="bottom") + 
  scale_color_manual(values=COL_GROUP_PAIRED, guide_legend(title = "")) +
  guides(fill=guide_legend(title="")) +
  labs(y="Expression Correlation [R]", x="")

ggsave(p, 
       file = paste0(
         outPrefix, 
        ".expCor_by_genotype_and_groups_sameTAD_Dixon_mESC.boxplot.pdf"), 
       w = 7, h = 7)


#-------------------------------------------------------------------------------
# All gene pairs: cor by Genotype, same TAD , TAD source, and GRB group 
#-------------------------------------------------------------------------------

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

# select only correlation across replicates
subDF <- subset(tidyDF, tidyDF$conditions == CONDITIONS)


p <- ggplot(subDF, aes(x = TAD, y = expCor, color = genotype)) + 
  geom_boxplot(lwd = 1.5) + 
  facet_grid(TADtype ~ study + tissue) + 
  theme_bw() + 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "bottom") + 
  scale_color_manual(values = COL_ALL_GENOTYPE, guide_legend(title = "")) +
  guides(fill = guide_legend(title = "")) +
  labs(y = "Expression Correlation [R]", x = "")

# geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=1.1), data=pvalDF, size=5) +
# geom_text(aes(label=paste0("n=",n), y=-1), data=nDF) +
# geom_text(aes(label=paste0("N=",nAll), y=-1.2), data=nDF)


ggsave(p, file = paste0(outPrefix, 
      ".all_groups.expCor_by_genotype_sameTAD_tadSource_GRB.boxplot.pdf"), 
      w = 7, h = 7)

#===============================================================================
# plot correlation in conditions against each other
#===============================================================================
# select expression data 
subDF <- subset(tidyDF, 
                singleTadSource &
                  tidyDF$cell == "Macrophages" & 
                  tidyDF$conditions == CONDITIONS)

# add two separate columns for WT and KO exp cor
corDF <- subDF %>%
  spread(key=genotype, value=expCor) 

# duplicate corDF to add group for all gene pairs 
nSub <- nrow(corDF)
allDF <- rbind(corDF, corDF)

levels(allDF$group) <- c(levels(allDF$group), "all")
allDF$group[seq(nSub + 1, 2 * nSub)] <- "all"

# assign names to colors
names(COL_GROUP) <- levels(allDF$group)

#-------------------------------------------------------------------------------
# correlation by gentoype and groups
#-------------------------------------------------------------------------------

# as density
p <- ggplot(allDF, aes(x = WT, y = Rad21KO)) +
  # geom_point(alpha = 0.1) + 
  #geom_density2d() + 
  # stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE) +
  # stat_density2d(geom = "raster", aes(alpha = ..density.., fill = group), contour = FALSE) +
  stat_density2d(geom = "tile", aes(fill = group, alpha = ..density..),
                 contour = FALSE) +
  scale_fill_manual(values = COL_GROUP) +
  facet_wrap(~ group) +
  #stat_density2d(aes(fill = ..level..), geom = "polygon") + 
  theme_bw() + 
  theme(text = element_text(size = 15), 
        axis.text.x=element_text(angle = 45, hjust = 1), 
        legend.position = "bottom") + 
  xlab("Correlation in WT [R]") + ylab("Correlation in Rad21 KO [R]")

ggsave(p, file=paste0(outPrefix, 
                      ".expCor_WTvsKO_by_group.pdf"), w=7, h=7)

#-------------
# as dotplot
#-------------
p <- ggplot(allDF, aes(x=WT, y=Rad21KO)) +
  geom_point(aes(color=group), alpha=0.25, size=.01) + 
  #geom_density2d() + 
  # stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE) + 
  # stat_density2d(geom = "raster", aes(alpha = ..density.., fill=group), contour = FALSE) +
  # stat_density2d(geom="tile", aes(fill = group, alpha=..density..), 
  #                contour=FALSE) + 
  scale_color_manual(values=COL_GROUP) +
  facet_wrap(~group) +
  #stat_density2d(aes(fill = ..level..), geom = "polygon") + 
  theme_bw() + 
  theme(text = element_text(size=15), 
        axis.text.x=element_text(angle = 45, hjust = 1), 
        legend.position="bottom") + 
  xlab("Correlation in WT [R]") + ylab("Correlation in Rad21 KO [R]")

ggsave(p, file = paste0(outPrefix, 
                      ".expCor_WTvsKO_by_group.dotplot.png"), w=7, h=7)

#-------------------------------------------------------------------------------
# correlation by gentoype and inTAD
#-------------------------------------------------------------------------------

p <- ggplot(corDF, aes(x=WT, y=Rad21KO)) +
  stat_density2d(geom="tile", aes(fill = TAD, alpha=..density..), 
                 contour=FALSE) + 
  scale_fill_manual(values=COL_TAD_GENOTYPE[c(1,3)]) +
  facet_grid(.~TAD) +
  theme_bw() + 
  theme(text = element_text(size=15), 
        axis.text.x=element_text(angle = 45, hjust = 1), 
        legend.position="bottom") + 
  xlab("Correlation in WT [R]") + 
  ylab("Correlation in Rad21 KO [R]")

ggsave(p, file=paste0(outPrefix, 
                      ".expCor_WTvsKO_by_TAD.pdf"), w=7, h=3.5)

#-------------
# as dotplot
#-------------

p <- ggplot(corDF, aes(x=WT, y=Rad21KO)) +
  geom_point(aes(color=TAD), alpha=0.25, size=.01) + 
  scale_color_manual(values=COL_TAD_GENOTYPE[c(1,3)]) +
  facet_grid(.~TAD) +
  theme_bw() + 
  theme(text = element_text(size=15), 
        axis.text.x=element_text(angle = 45, hjust = 1), 
        legend.position="bottom") + 
  xlab("Correlation in WT [R]") + 
  ylab("Correlation in Rad21 KO [R]")

ggsave(p, file=paste0(outPrefix, 
                      ".expCor_WTvsKO_by_TAD.png"), w=7, h=3.5)



#===============================================================================
# plot change of expression values upon knock-down
#===============================================================================
# select expression data 
subDF <- subset(tidyDF, 
                singleTadSource &
                tidyDF$cell == "Macrophages" & 
                  tidyDF$conditions == CONDITIONS)

# compute difference in expression between KO and WT
deltaDF <- subDF %>%
  spread(key=genotype, value=expCor) %>%
  mutate(deltaExpCor = Rad21KO - WT) 

#-------------------------------------------------------------------------------
# plot corelation change for all groups
#-------------------------------------------------------------------------------

groupDF <- deltaDF %>%
  group_by(TAD) %>%
  summarise(
    mean=mean(deltaExpCor, na.rm=TRUE),
    median=median(deltaExpCor, na.rm=TRUE),
    sd=sd(deltaExpCor, na.rm=TRUE),
    N=n(),
    n=sum(!is.na(deltaExpCor))
  )

pVal <- wilcox.test(deltaExpCor ~ TAD, data=deltaDF)$p.value

p <- ggplot(deltaDF, aes(x=TAD, y=deltaExpCor, color=TAD)) + 
  geom_boxplot(lwd=1.5) + 
  facet_grid(TADtype ~ study+tissue) +
  geom_label(aes(y=-2, label=paste0(
    "avg.=", signif(mean, 2), "\n",
    "med.=", signif(median, 2), "\n",
    "n=", n
  )), data=groupDF, vjust = "bottom", hjust = "center") +
  geom_text(aes(label=paste0("p=", signif(pVal, 2)), x=1.5, y=1.5), color="black", data=data.frame()) +
  theme_bw() + 
  theme(text = element_text(size=15), 
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position="none") + 
  scale_color_manual(values=COL_TAD_GENOTYPE[c(1,3)]) +
  guides(fill=guide_legend(title="")) +
  labs(y="Expression Correlation Change [Rad21KO - WT]", x="")


ggsave(p, file=paste0(outPrefix, 
                      ".DeltaExpCor_by_sameTAD_Dixon_mESC.boxplot.pdf"), w=3.5, h=7)

#-------------------------------------------------------------------------------
# plot corelation change by groups
#-------------------------------------------------------------------------------
groupDF <- deltaDF %>%
  group_by(TAD, group) %>%
  summarise(
    mean=mean(deltaExpCor, na.rm=TRUE),
    median=median(deltaExpCor, na.rm=TRUE),
    sd=sd(deltaExpCor, na.rm=TRUE),
    N=n(),
    n=sum(!is.na(deltaExpCor))
  )


# calculate p-values for differnece in delta by TAD
pvalDF <- deltaDF %>%
  group_by(group) %>%
  do(w = wilcox.test(deltaExpCor ~ TAD, data=.)) %>%
  summarise(group, p=w$p.value)

COL=mixColors(COL_GROUP, COL_TAD_GENOTYPE[c(1,3)], 0.4)

p <- ggplot(deltaDF, aes(x=TAD, y=deltaExpCor, color=interaction(TAD,group))) + 
  geom_boxplot(lwd=1.5) + 
  facet_grid( ~ group) +
  geom_label(aes(y=-2, label=paste0(
    "avg.=", signif(mean, 2), "\n",
    "med.=", signif(median, 2), "\n",
    "n=", n
  )), data=groupDF, vjust = "bottom", hjust = "center") +
  geom_text(aes(label=paste0("p=", signif(p, 2)), x=1.5, y=1.5), color="black", data=pvalDF) +
  theme_bw() + 
  theme(text = element_text(size=15), 
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position="none") + 
  scale_color_manual(values=COL) +
  guides(fill=guide_legend(title="")) +
  labs(y="Expression Correlation Change [Rad21KO - WT]", x="")


ggsave(p, file=paste0(outPrefix, 
                      ".DeltaExpCor_by_sameTAD_and_group_Dixon_mESC.boxplot.pdf"), w=7, h=7)



#-------------------------------------------------------------------------------
# plot corelation change for by tad source
#-------------------------------------------------------------------------------
# select expression data 
subDF <- subset(tidyDF, 
                  tidyDF$cell == "Macrophages" & 
                  tidyDF$conditions == CONDITIONS)

# compute difference in expression between KO and WT
deltaDF <- subDF %>%
  spread(key=genotype, value=expCor) %>%
  mutate(deltaExpCor = Rad21KO - WT)


groupDF <- deltaDF %>%
  group_by(TAD, study, tissue, TADtype) %>%
  summarise(
    mean=mean(deltaExpCor, na.rm=TRUE),
    median=median(deltaExpCor, na.rm=TRUE),
    sd=sd(deltaExpCor, na.rm=TRUE),
    N=n(),
    n=sum(!is.na(deltaExpCor))
  )

# calculate p-values for differnece in delta by TAD
pvalDF <- deltaDF %>%
  group_by(study, tissue, TADtype) %>%
  do(w = wilcox.test(deltaExpCor ~ TAD, data=.)) %>%
  summarise(study, tissue, TADtype, p=w$p.value)



p <- ggplot(deltaDF, aes(x=TAD, y=deltaExpCor, col=TAD)) + 
  geom_boxplot(lwd=1) + 
  facet_grid(TADtype ~ study+tissue) + 
  geom_label(aes(y=-2, label=paste0(
    "avg.=", signif(mean, 2), "\n",
    "med.=", signif(median, 2), "\n",
    "n=", n
  )), data=groupDF, vjust = "bottom", hjust = "center", size=2) +
  geom_text(aes(label=paste0("p=", signif(p, 2)), x=1.5, y=1.5), color="black", data=pvalDF) +
  theme_bw() + 
  theme(text = element_text(size=15), 
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position="none") + 
  scale_color_manual(values=COL_TAD_GENOTYPE[c(1,3)]) +
  guides(fill=guide_legend(title="")) +
  labs(y="Expression Correlation Change [Rad21KO - WT]", x="")


ggsave(p, file=paste0(outPrefix, 
                      ".DeltaExpCor_by_sameTAD_group_tadSource_GRB.boxplot.pdf"), w=7, h=7)


#===============================================================================
# plot expression correlation by distance
#===============================================================================

# select only correlation across replicates
subDF <- subset(tidyDF, singleTadSource & 
                  tidyDF$conditions == CONDITIONS)

ALPHA=.3


p <- ggplot(subDF, aes(x=dist, y=expCor, linetype=genotype, fill=TAD, col=TAD)) + 
  # geom_point(alpha=ALPHA, size=.5, shape=20)  +
  geom_smooth(alpha=ALPHA) + 
  guides(col=guide_legend(title=""), fill=guide_legend(title=""), 
         linetype=guide_legend(title="")) + 
  scale_color_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + 
  scale_fill_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + 
  theme_bw() + 
  theme(legend.position="bottom") +
  labs(x="Distance between gene pairs [kb]", y="Expression Correlation [R]")

    # + guides(linetype=FALSE)
  
ggsave(p, 
       file=paste0(outPrefix, 
        ".all_groups.expCor_by_distance_genotype_sameTAD_Dixon_mESC.boxplot.pdf"), 
       w=7, h=7)


p <- ggplot(subDF, aes(x=dist, y=expCor, linetype=genotype, fill=Boundary, col=Boundary)) + 
  # geom_point(alpha=ALPHA, size=.5, shape=20)  +
  geom_smooth(alpha=ALPHA) + 
  guides(col=guide_legend(title=""), fill=guide_legend(title=""),
         linetype=guide_legend(title="")) + 
  scale_color_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + 
  scale_fill_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + 
  theme_bw() + 
  theme(legend.position="bottom") +
  labs(x="Distance between gene pairs [kb]", y="Expression Correlation [R]")

ggsave(p, 
       file=paste0(outPrefix, 
                   ".all_groups.expCor_by_distance_genotype_NotSeparated_Dixon_mESC.boxplot.pdf"), 
       w=7, h=7)

#-------------------------------------------------------------------------------
# Expression by distance for all TAD sources
#-------------------------------------------------------------------------------
# select only correlation across replicates
subDF <- subset(tidyDF, tidyDF$conditions == CONDITIONS)


ALPHA=.3


p <- ggplot(subDF, aes(x=dist, y=expCor, linetype=genotype, fill=TAD, col=TAD)) + 
  # geom_point(alpha=ALPHA, size=.5, shape=20)  +
  geom_smooth(alpha=ALPHA) + 
  facet_grid(TADtype ~ study+tissue) + 
  guides(col=guide_legend(title=""), fill=guide_legend(title="")) + 
  scale_color_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + scale_fill_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + 
  theme_bw() + 
  theme(legend.position="bottom") +
  labs(x="Distance between gene pairs [kb]", y="Expression Correlation [R]")

# + guides(linetype=FALSE)

ggsave(p, 
       file=paste0(outPrefix, 
                   ".all_groups.expCor_by_distance_genotype_sameTAD_tadSource_GRB.boxplot.pdf"), 
       w=7, h=7)


p <- ggplot(subDF, aes(x=dist, y=expCor, linetype=genotype, fill=Boundary, col=Boundary)) + 
  # geom_point(alpha=ALPHA, size=.5, shape=20)  +
  geom_smooth(alpha=ALPHA) + 
  facet_grid(TADtype ~ study+tissue) + 
  guides(col=guide_legend(title=""), fill=guide_legend(title="")) + 
  scale_color_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + scale_fill_manual(values = COL_TAD_GENOTYPE[c(1,3)]) + 
  theme_bw() + 
  theme(legend.position="bottom") +
  labs(x="Distance between gene pairs [kb]", y="Expression Correlation [R]")

ggsave(p, 
       file=paste0(outPrefix, 
                   ".all_groups.expCor_by_distance_genotype_NotSeparated_tadSource_GRB.boxplot.pdf.pdf"), 
       w=7, h=7)

