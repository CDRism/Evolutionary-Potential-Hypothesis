#Effects of testosterone on gene expression are concordant between sexes but divergent across species of Sceloporus lizards

library(edgeR)
library(DESeq2)
library(plyr)
library(dplyr)
library(ggplot2)
library(rospca)
library(cowplot)

#Used throughout for consistency
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plot_theme = theme_bw()+ 
  theme(axis.text=element_text(size=20, color = "black"), axis.title=element_text(size=20,face="bold"))+
  theme(legend.text = element_blank(), legend.title = element_blank())+
  theme(axis.ticks.length = unit(0.4, "cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.line = element_line(size = 1.5))+
  theme(legend.position="none")+
  theme(panel.grid.major=element_line(colour="white"), panel.grid.minor = element_line(colour = "white"))+
  theme(panel.border=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0),"cm"))


####Figure 1 and 2: Sex-specific and shared responses to testosterone####
data = read.csv("ReadCounts.csv")[,-c(3:6)]
#to remove the X chromosome:
chrs = c("NC_056522.1","NC_056523.1","NC_056524.1","NC_056525.1","NC_056526.1","NC_056527.1","NC_056528.1","NC_056529.1","NC_056530.1","NC_056532.1")
data$Chr = substr(data$Chr,1,11) #keep only the chromosome information
data = data[data$Chr %in% chrs,] #keep only the chromosomes of interest
data = data[complete.cases(data),]
rownames(data) = data[,1]
data = data[,-c(1:2)]
data = data[,-c(4,22)] #removes two bad merriami

####undulatus sex specificity####
und.data = data[,c(22:46)]

#get undulatus group data
spec.data = read.csv("morphological_data.csv")
spec.data = spec.data[spec.data$ID %in% colnames(data),]
group = paste(spec.data$Species, spec.data$Sex, spec.data$Treatment, sep = "")
group.und = group[22:46]
treat.und = spec.data[22:46,4]
sex.und = spec.data[22:46,3]

#create models
DEGs = DGEList(counts = und.data, group = group.und)
design =model.matrix(~0+group.und, data = DEGs$samples)

keep = filterByExpr(DEGs, design = design, group = group.und)
DEGs = DEGs[keep,] #keeps 13891 genes
und.data = und.data[rownames(und.data) %in% rownames(DEGs),]
DEGs = calcNormFactors(DEGs)
logcpm = cpm(DEGs, log = T)
DEGs = estimateDisp(DEGs, design)
DEGs$common.dispersion #0.172
fit_robust_auto = glmQLFit(DEGs, design, robust = T)

#males only
contrasts_auto = makeContrasts(und.male.TvC = group.undSCUNMTEST - group.undSCUNMCONT, levels = design)
QLF.und.male.TvC = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.male.TvC"])
de = decideTestsDGE(QLF.und.male.TvC, p = 0.05)
summary(de) #22 to testosterone, 19 to controls

logfc.und.male = data.frame(topTags(QLF.und.male.TvC, n = Inf))
logfc.TEST.und.male = subset(logfc.und.male, logFC > 0 & FDR < 0.05)
logfc.CONT.und.male = subset(logfc.und.male, logFC < 0 & FDR < 0.05)

und.male = ggplot(logfc.und.male, aes(x = logFC, y = -log10(PValue)))+
  geom_point(color = "lightgray", size = 3)+
  geom_point(fill = ifelse(rownames(logfc.und.male) %in% rownames(logfc.TEST.und.male), "dodgerblue3",
                           ifelse(rownames(logfc.und.male) %in% rownames(logfc.CONT.und.male), "indianred1", NA)),
             size = 3,
             alpha = ifelse(rownames(logfc.und.male) %in% rownames(logfc.TEST.und.male) |
                              rownames(logfc.und.male) %in% rownames(logfc.CONT.und.male), 1, 0),
             pch = 21)+
  scale_x_continuous(breaks = c(-11,-5.5,0,5.5,11))+
  coord_cartesian(ylim = c(0,12), xlim = c(-11,11), expand = F)+
  xlab(expression(bold("log"["2"]*"FC")))+
  ylab(expression(bold("-log"["10"]*"("*bolditalic(P)*"-Value)")))+
  plot_theme+theme(plot.margin = unit(c(0.25,0.1,0.1,0.1),"cm"))

#females only
contrasts_auto = makeContrasts(und.female.TvC = group.undSCUNFTEST - group.undSCUNFCONT, levels = design)
QLF.und.female.TvC = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.female.TvC"])
de = decideTestsDGE(QLF.und.female.TvC, p = 0.05)
summary(de) #19 to testosterone, 6 to controls

logfc.und.female = data.frame(topTags(QLF.und.female.TvC, n = Inf))
logfc.TEST.und.female = subset(logfc.und.female, logFC > 0 & FDR < 0.05)
logfc.CONT.und.female = subset(logfc.und.female, logFC < 0 & FDR < 0.05)

und.female = ggplot(logfc.und.female, aes(x = logFC, y = -log10(PValue)))+
  geom_point(color = "lightgray", size = 3)+
  geom_point(fill = ifelse(rownames(logfc.und.female) %in% rownames(logfc.TEST.und.female), "dodgerblue3",
                           ifelse(rownames(logfc.und.female) %in% rownames(logfc.CONT.und.female), "indianred1", NA)),
             size = 3,
             alpha = ifelse(rownames(logfc.und.female) %in% rownames(logfc.TEST.und.female) |
                              rownames(logfc.und.female) %in% rownames(logfc.CONT.und.female), 1, 0),
             pch = 21)+
  scale_x_continuous(breaks = c(-11,-5.5,0,5.5,11))+
  coord_cartesian(ylim = c(0,12), xlim = c(-11,11), expand = F)+
  xlab(expression(bold("log"["2"]*"FC")))+
  ylab(expression(bold("-log"["10"]*"("*bolditalic(P)*"-Value)")))+
  plot_theme+theme(plot.margin = unit(c(0.25,0.1,0.1,0.1),"cm"))

#does number of DEGs differ?
#overall
expected = mean(c((22+19),(19+6)))
o1 = (22 + 19) - expected
o2 = (19 + 6) - expected
o1_standardized = o1^2/expected
o2_standardized = o2^2/expected
chi = o1_standardized + o2_standardized
chi > 3.841

#testosterone
expected = mean(c(22,19))
o1 = 22 - expected
o2 = 19 - expected
o1_standardized = o1^2/expected
o2_standardized = o2^2/expected
chi = o1_standardized + o2_standardized
chi > 3.841

#control
expected = mean(c(19,6))
o1 = 19 - expected
o2 = 6 - expected
o1_standardized = o1^2/expected
o2_standardized = o2^2/expected
chi = o1_standardized + o2_standardized
chi > 3.841

#both sexes
contrasts_auto = makeContrasts(und.TvC = (group.undSCUNFTEST+group.undSCUNMTEST)/2 - (group.undSCUNFCONT+group.undSCUNMCONT)/2, levels = design)
QLF.und.TvC = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.TvC"])
de = decideTestsDGE(QLF.und.TvC, p = 0.05)
summary(de) #126 to testosterone, 157 to controls

logfc.und = data.frame(topTags(QLF.und.TvC, n = Inf))
logfc.TEST.und = subset(logfc.und, logFC > 0 & FDR < 0.05)
logfc.CONT.und = subset(logfc.und, logFC < 0 & FDR < 0.05)

und.both = ggplot(logfc.und, aes(x = logFC, y = -log10(PValue)))+
  geom_point(color = "lightgray", size = 3)+
  geom_point(fill = ifelse(rownames(logfc.und) %in% rownames(logfc.TEST.und), "dodgerblue3",
                           ifelse(rownames(logfc.und) %in% rownames(logfc.CONT.und), "indianred1", NA)),
             size = 3,
             alpha = ifelse(rownames(logfc.und) %in% rownames(logfc.TEST.und) |
                              rownames(logfc.und) %in% rownames(logfc.CONT.und), 1, 0),
             pch = 21)+
  scale_x_continuous(breaks = c(-11,-5.5,0,5.5,11))+
  scale_y_continuous(breaks = seq(0,16,by = 4))+
  coord_cartesian(ylim = c(0,16), xlim = c(-11,11), expand = F)+
  xlab(expression(bold("log"["2"]*"FC")))+
  ylab(expression(bold("-log"["10"]*"("*bolditalic(P)*"-Value)")))+
  plot_theme+theme(plot.margin = unit(c(0.25,0.1,0.1,0.1),"cm"))

#interaction
contrasts_auto = makeContrasts(und.int = (group.undSCUNMTEST - group.undSCUNMCONT) - (group.undSCUNFTEST - group.undSCUNFCONT), levels = design)
QLF.und.int = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.int"])
de = decideTestsDGE(QLF.und.int, p = 0.05)
summary(de) #none

logfc.und.int = data.frame(topTags(QLF.und.int, n = Inf))
logfc.und.up = subset(logfc.und.int, logFC > 0 & FDR < 0.05)
logfc.und.down = subset(logfc.und.int, logFC < 0 & FDR < 0.05)

data.merged = merge(logfc.und.male, logfc.und.female, by = 0) #0 is rownames
data.merged$gene = data.merged$Row.names

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$Row.names))
colnames(fig.data) = c("Male", "Female", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$Male = as.numeric(fig.data$Male)
fig.data$Female = as.numeric(fig.data$Female)
fig.data$density = get_density(fig.data$Male, fig.data$Female, n = 1000)

und.int.plot = ggplot(fig.data, aes(x = Male, y = Female))+
  geom_point(aes(x = Male, y = Female, col = density))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.TEST.und),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.CONT.und),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.und.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.und.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  xlab(expression(bolditalic("S. undulatus")*bold(" Male log"["2"]*"FC")))+
  ylab(expression(bolditalic("S. undulatus ")*bold(" Female log"["2"]*"FC")))+
  plot_theme
summary(lm(fig.data$Female ~ fig.data$Male)) #r2 = 0.212
sqrt(0.2122) #0.460

#controls only
contrasts_auto = makeContrasts(und.cont = group.undSCUNMCONT - group.undSCUNFCONT, levels = design)
QLF.und.cont = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.cont"])
de = decideTestsDGE(QLF.und.cont, p = 0.05)
summary(de) #1 higher in males

logfc.und.cont = data.frame(topTags(QLF.und.cont, n = Inf))
logfc.CONT.und.sexes = subset(logfc.und.cont, logFC > 0 & FDR < 0.05)

#test only
contrasts_auto = makeContrasts(und.test = group.undSCUNMTEST - group.undSCUNFTEST, levels = design)
QLF.und.test = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.test"])
de = decideTestsDGE(QLF.und.test, p = 0.05)
summary(de) #no differences


####virgatus sex specificity####
virg.data = data[,c(47:70)]
spec.data = read.csv("morphological_data.csv")
spec.data = spec.data[spec.data$ID %in% colnames(data),]
group = paste(spec.data$Species, spec.data$Sex, spec.data$Treatment, sep = "")
group.virg = group[c(47:70)]

DEGs = DGEList(counts = virg.data, group = group.virg)
design = model.matrix(~0+group.virg, data = DEGs$samples)
keep = filterByExpr(DEGs, design = design, group = group.virg)
DEGs = DEGs[keep,] #13,772 genes
virg.data = virg.data[rownames(virg.data) %in% rownames(DEGs),]

DEGs = calcNormFactors(DEGs)
logcpm = cpm(DEGs, log = T)
DEGs = estimateDisp(DEGs, design)
DEGs$common.dispersion #0.181
fit_robust_auto = glmQLFit(DEGs, design, robust = T)

#males only
contrasts_auto = makeContrasts(virg.male.TvC = group.virgSCVIMTEST - group.virgSCVIMCONT, levels = design)
QLF.virg.male.TvC = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.male.TvC"])
de = decideTestsDGE(QLF.virg.male.TvC, p = 0.05)
summary(de) #69 to testosterone, 12 to controls

logfc.virg.male = data.frame(topTags(QLF.virg.male.TvC, n = Inf))
logfc.TEST.virg.male = subset(logfc.virg.male, logFC > 0 & FDR < 0.05)
logfc.CONT.virg.male = subset(logfc.virg.male, logFC < 0 & FDR < 0.05)

virg.male = ggplot(logfc.virg.male, aes(x = logFC, y = -log10(PValue)))+
  geom_point(color = "lightgray", size = 3)+
  geom_point(fill = ifelse(rownames(logfc.virg.male) %in% rownames(logfc.TEST.virg.male), "dodgerblue3",
                           ifelse(rownames(logfc.virg.male) %in% rownames(logfc.CONT.virg.male), "indianred1", NA)),
             size = 3,
             alpha = ifelse(rownames(logfc.virg.male) %in% rownames(logfc.TEST.virg.male) |
                              rownames(logfc.virg.male) %in% rownames(logfc.CONT.virg.male), 1, 0),
             pch = 21)+
  scale_x_continuous(breaks = c(-11,-5.5,0,5.5,11))+
  scale_y_continuous(breaks = seq(0,12,by=3))+
  coord_cartesian(ylim = c(0,12), xlim = c(-11,11), expand = F)+
  xlab(expression(bold("log"["2"]*"FC")))+
  ylab(expression(bold("-log"["10"]*"("*bolditalic(P)*"-Value)")))+
  plot_theme+theme(plot.margin = unit(c(0.25,0.1,0.1,0.1),"cm"))

#females only
contrasts_auto = makeContrasts(virg.female.TvC = group.virgSCVIFTEST - group.virgSCVIFCONT, levels = design)
QLF.virg.female.TvC = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.female.TvC"])
de = decideTestsDGE(QLF.virg.female.TvC, p = 0.05)
summary(de) #27 to testosterone, 0 to controls

logfc.virg.female = data.frame(topTags(QLF.virg.female.TvC, n = Inf))
logfc.TEST.virg.female = subset(logfc.virg.female, logFC > 0 & FDR < 0.05)
logfc.CONT.virg.female = subset(logfc.virg.female, logFC < 0 & FDR < 0.05)

virg.female = ggplot(logfc.virg.female, aes(x = logFC, y = -log10(PValue)))+
  geom_point(color = "lightgray", size = 3)+
  geom_point(fill = ifelse(rownames(logfc.virg.female) %in% rownames(logfc.TEST.virg.female), "dodgerblue3",
                           ifelse(rownames(logfc.virg.female) %in% rownames(logfc.CONT.virg.female), "indianred1", NA)),
             size = 3,
             alpha = ifelse(rownames(logfc.virg.female) %in% rownames(logfc.TEST.virg.female) |
                              rownames(logfc.virg.female) %in% rownames(logfc.CONT.virg.female), 1, 0),
             pch = 21)+
  scale_x_continuous(breaks = c(-11,-5.5,0,5.5,11))+
  scale_y_continuous(breaks = seq(0,12,by=3))+
  coord_cartesian(ylim = c(0,12), xlim = c(-11,11), expand = F)+
  xlab(expression(bold("log"["2"]*"FC")))+
  ylab(expression(bold("-log"["10"]*"("*bolditalic(P)*"-Value)")))+
  plot_theme+theme(plot.margin = unit(c(0.25,0.1,0.1,0.1),"cm"))

#does number of DEGs differ?
#overall
expected = mean(c((69+12),(27+0)))
o1 = (69 + 12) - expected
o2 = (27 + 0) - expected
o1_standardized = o1^2/expected
o2_standardized = o2^2/expected
chi = o1_standardized + o2_standardized
chi > 3.841

#testosterone
expected = mean(c(69,27))
o1 = 69 - expected
o2 = 27 - expected
o1_standardized = o1^2/expected
o2_standardized = o2^2/expected
chi = o1_standardized + o2_standardized
chi > 3.841

#control
expected = mean(c(12,0))
o1 = 12 - expected
o2 = 0 - expected
o1_standardized = o1^2/expected
o2_standardized = o2^2/expected
chi = o1_standardized + o2_standardized
chi > 3.841

#both sexes
contrasts_auto = makeContrasts(virg.TvC = (group.virgSCVIFTEST+group.virgSCVIMTEST)/2 - (group.virgSCVIFCONT+group.virgSCVIMCONT)/2, levels = design)
QLF.virg.TvC = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.TvC"])
de = decideTestsDGE(QLF.virg.TvC, p = 0.05)
summary(de) #465 to testosterone, 159 to controls

logfc.virg = data.frame(topTags(QLF.virg.TvC, n = Inf))
logfc.TEST.virg = subset(logfc.virg, logFC > 0 & FDR < 0.05)
logfc.CONT.virg = subset(logfc.virg, logFC < 0 & FDR < 0.05)

virg.both = ggplot(logfc.virg, aes(x = logFC, y = -log10(PValue)))+
  geom_point(color = "lightgray", size = 3)+
  geom_point(fill = ifelse(rownames(logfc.virg) %in% rownames(logfc.TEST.virg), "dodgerblue3",
                           ifelse(rownames(logfc.virg) %in% rownames(logfc.CONT.virg), "indianred1", NA)),
             size = 3,
             alpha = ifelse(rownames(logfc.virg) %in% rownames(logfc.TEST.virg) |
                              rownames(logfc.virg) %in% rownames(logfc.CONT.virg), 1, 0),
             pch = 21)+
  scale_x_continuous(breaks = c(-11,-5.5,0,5.5,11))+
  coord_cartesian(ylim = c(0,12), xlim = c(-11,11), expand = F)+
  xlab(expression(bold("log"["2"]*"FC")))+
  ylab(expression(bold("-log"["10"]*"("*bolditalic(P)*"-Value)")))+
  plot_theme+theme(plot.margin = unit(c(0.25,0.1,0.1,0.1),"cm"))

#interaction
contrasts_auto = makeContrasts(virg.int = (group.virgSCVIMTEST - group.virgSCVIMCONT) - (group.virgSCVIFTEST - group.virgSCVIFCONT), levels = design)
QLF.virg.int = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.int"])
de = decideTestsDGE(QLF.virg.int, p = 0.05)
summary(de) #none

logfc.virg.int = data.frame(topTags(QLF.virg.int, n = Inf))
logfc.virg.up = subset(logfc.virg.int, logFC > 0 & FDR < 0.05)
logfc.virg.down = subset(logfc.virg.int, logFC < 0 & FDR < 0.05)

data.merged = merge(logfc.virg.male, logfc.virg.female, by = 0) #0 is rownames
data.merged$gene = data.merged$Row.names

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$Row.names))
colnames(fig.data) = c("Male", "Female", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$Male = as.numeric(fig.data$Male)
fig.data$Female = as.numeric(fig.data$Female)
fig.data$density = get_density(fig.data$Male, fig.data$Female, n = 1000)

virg.int.plot = ggplot(fig.data, aes(x = Male, y = Female))+
  geom_point(aes(x = Male, y = Female, col = density))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.TEST.virg),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.CONT.virg),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virg.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virg.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$Female ~ fig.data$Male)) #r2 = 0.202
sqrt(0.2024) #0.450

#controls only
contrasts_auto = makeContrasts(virg.cont = group.virgSCVIMCONT - group.virgSCVIFCONT, levels = design)
QLF.virg.cont = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.cont"])
de = decideTestsDGE(QLF.virg.cont, p = 0.05)
summary(de) #no differences

#test only
contrasts_auto = makeContrasts(virg.test = group.virgSCVIMTEST - group.virgSCVIFTEST, levels = design)
QLF.virg.test = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.test"])
de = decideTestsDGE(QLF.virg.test, p = 0.05)
summary(de) #no differences


####merriami sex specificity####
mer.data = data[,c(1:21)]
spec.data = read.csv("morphological_data.csv")
spec.data = spec.data[spec.data$ID %in% colnames(data),]
group = paste(spec.data$Species, spec.data$Sex, spec.data$Treatment, sep = "")
group.mer = group[1:21]

DEGs = DGEList(counts = mer.data, group = group.mer)
design = model.matrix(~0+group.mer, data = DEGs$samples)
keep = filterByExpr(DEGs, design = design, group = group.mer)
DEGs = DEGs[keep,] #keeps 13036 genes
mer.data = data[rownames(mer.data) %in% rownames(DEGs),]
DEGs = calcNormFactors(DEGs)
logcpm = cpm(DEGs, log = T)
DEGs = estimateDisp(DEGs, design)
DEGs$common.dispersion #0.164
fit_robust_auto = glmQLFit(DEGs, design, robust = T)

#males only
contrasts_auto = makeContrasts(mer.male.TvC = group.merSCMEMTEST - group.merSCMEMCONT, levels = design)
QLF.mer.male.TvC = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.male.TvC"])
de = decideTestsDGE(QLF.mer.male.TvC, p = 0.05)
summary(de) #192 to testosterone, 156 to controls

logfc.mer.male = data.frame(topTags(QLF.mer.male.TvC, n = Inf))
logfc.TEST.mer.male = subset(logfc.mer.male, logFC > 0 & FDR < 0.05)
logfc.CONT.mer.male = subset(logfc.mer.male, logFC < 0 & FDR < 0.05)

mer.male = ggplot(logfc.mer.male, aes(x = logFC, y = -log10(PValue)))+
  geom_point(color = "lightgray", size = 3)+
  geom_point(fill = ifelse(rownames(logfc.mer.male) %in% rownames(logfc.TEST.mer.male), "dodgerblue3",
                           ifelse(rownames(logfc.mer.male) %in% rownames(logfc.CONT.mer.male), "indianred1", NA)),
             size = 3,
             alpha = ifelse(rownames(logfc.mer.male) %in% rownames(logfc.TEST.mer.male) |
                              rownames(logfc.mer.male) %in% rownames(logfc.CONT.mer.male), 1, 0),
             pch = 21)+
  scale_x_continuous(breaks = c(-11,-5.5,0,5.5,11))+
  scale_y_continuous(breaks = seq(0,12,by=3))+
  coord_cartesian(ylim = c(0,12), xlim = c(-11,11), expand = F)+
  xlab(expression(bold("log"["2"]*"FC")))+
  ylab(expression(bold("-log"["10"]*"("*bolditalic(P)*"-Value)")))+
  plot_theme+theme(plot.margin = unit(c(0.25,0.1,0.1,0.1),"cm"))

#females only
contrasts_auto = makeContrasts(mer.female.TvC = group.merSCMEFTEST - group.merSCMEFCONT, levels = design)
QLF.mer.female.TvC = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.female.TvC"])
de = decideTestsDGE(QLF.mer.female.TvC, p = 0.05)
summary(de) #79 to testosterone, 18 to controls

logfc.mer.female = data.frame(topTags(QLF.mer.female.TvC, n = Inf))
logfc.TEST.mer.female = subset(logfc.mer.female, logFC > 0 & FDR < 0.05)
logfc.CONT.mer.female = subset(logfc.mer.female, logFC < 0 & FDR < 0.05)

mer.female = ggplot(logfc.mer.female, aes(x = logFC, y = -log10(PValue)))+
  geom_point(color = "lightgray", size = 3)+
  geom_point(fill = ifelse(rownames(logfc.mer.female) %in% rownames(logfc.TEST.mer.female), "dodgerblue3",
                           ifelse(rownames(logfc.mer.female) %in% rownames(logfc.CONT.mer.female), "indianred1", NA)),
             size = 3,
             alpha = ifelse(rownames(logfc.mer.female) %in% rownames(logfc.TEST.mer.female) |
                              rownames(logfc.mer.female) %in% rownames(logfc.CONT.mer.female), 1, 0),
             pch = 21)+
  scale_x_continuous(breaks = c(-11,-5.5,0,5.5,11))+
  scale_y_continuous(breaks = seq(0,12,by=3))+
  coord_cartesian(ylim = c(0,12), xlim = c(-11,11), expand = F)+
  xlab(expression(bold("log"["2"]*"FC")))+
  ylab(expression(bold("-log"["10"]*"("*bolditalic(P)*"-Value)")))+
  plot_theme+theme(plot.margin = unit(c(0.25,0.1,0.1,0.1),"cm"))

#does number of DEGs differ?
#overall
expected = mean(c((156+192),(18+79)))
o1 = (156 + 192) - expected
o2 = (18 + 79) - expected
o1_standardized = o1^2/expected
o2_standardized = o2^2/expected
chi = o1_standardized + o2_standardized
chi > 3.841

#testosterone
expected = mean(c(192,79))
o1 = 192 - expected
o2 = 79 - expected
o1_standardized = o1^2/expected
o2_standardized = o2^2/expected
chi = o1_standardized + o2_standardized
chi > 3.841

#control
expected = mean(c(156,18))
o1 = 156 - expected
o2 = 18 - expected
o1_standardized = o1^2/expected
o2_standardized = o2^2/expected
chi = o1_standardized + o2_standardized
chi > 3.841

#both sexes
contrasts_auto = makeContrasts(mer.TvC = (group.merSCMEFTEST+group.merSCMEMTEST)/2 - (group.merSCMEFCONT+group.merSCMEMCONT)/2, levels = design)
QLF.mer.TvC = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.TvC"])
de = decideTestsDGE(QLF.mer.TvC, p = 0.05)
summary(de) #467 to testosterone, 393 to controls

logfc.mer = data.frame(topTags(QLF.mer.TvC, n = Inf))
logfc.TEST.mer = subset(logfc.mer, logFC > 0 & FDR < 0.05)
logfc.CONT.mer = subset(logfc.mer, logFC < 0 & FDR < 0.05)

mer.both = ggplot(logfc.mer, aes(x = logFC, y = -log10(PValue)))+
  geom_point(color = "lightgray", size = 3)+
  geom_point(fill = ifelse(rownames(logfc.mer) %in% rownames(logfc.TEST.mer), "dodgerblue3",
                           ifelse(rownames(logfc.mer) %in% rownames(logfc.CONT.mer), "indianred1", NA)),
             size = 3,
             alpha = ifelse(rownames(logfc.mer) %in% rownames(logfc.TEST.mer) |
                              rownames(logfc.mer) %in% rownames(logfc.CONT.mer), 1, 0),
             pch = 21)+
  scale_x_continuous(breaks = c(-11,-5.5,0,5.5,11))+
  scale_y_continuous(breaks = seq(0,12,by=3))+
  coord_cartesian(ylim = c(0,12), xlim = c(-11,11), expand = F)+
  xlab(expression(bold("log"["2"]*"FC")))+
  ylab(expression(bold("-log"["10"]*"("*bolditalic(P)*"-Value)")))+
  plot_theme+theme(plot.margin = unit(c(0.25,0.1,0.1,0.1),"cm"))

#interaction
contrasts_auto = makeContrasts(mer.int = (group.merSCMEMTEST - group.merSCMEMCONT) - (group.merSCMEFTEST - group.merSCMEFCONT), levels = design)
QLF.mer.int = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.int"])
de = decideTestsDGE(QLF.mer.int, p = 0.05)
summary(de) #none

logfc.mer.int = data.frame(topTags(QLF.mer.int, n = Inf))
logfc.mer.up = subset(logfc.mer.int, logFC > 0 & FDR < 0.05)
logfc.mer.down = subset(logfc.mer.int, logFC < 0 & FDR < 0.05)

data.merged = merge(logfc.mer.male, logfc.mer.female, by = 0) #0 is rownames
data.merged$gene = data.merged$Row.names

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$Row.names))
colnames(fig.data) = c("Male", "Female", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$Male = as.numeric(fig.data$Male)
fig.data$Female = as.numeric(fig.data$Female)
fig.data$density = get_density(fig.data$Male, fig.data$Female, n = 1000)

mer.int.plot = ggplot(fig.data, aes(x = Male, y = Female))+
  geom_point(aes(x = Male, y = Female, col = density))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.TEST.mer),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.CONT.mer),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.mer.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.mer.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$Female ~ fig.data$Male)) #r2 = 0.304
sqrt(0.3038) #0.551

#controls only
contrasts_auto = makeContrasts(mer.cont = group.merSCMEMCONT - group.merSCMEFCONT, levels = design)
QLF.mer.cont = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.cont"])
de = decideTestsDGE(QLF.mer.cont, p = 0.05)
summary(de) #no differences

#test only
contrasts_auto = makeContrasts(mer.test = group.merSCMEMTEST - group.merSCMEFTEST, levels = design)
QLF.mer.test = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.test"])
de = decideTestsDGE(QLF.mer.test, p = 0.05)
summary(de) #no differences


####Figure 3: cross-species interactions####
#note - this requires an omnibus model to do correctly.
data = read.csv("ReadCounts.csv")[,-c(3:6)]
chrs = c("NC_056522.1","NC_056523.1","NC_056524.1","NC_056525.1","NC_056526.1","NC_056527.1","NC_056528.1","NC_056529.1","NC_056530.1","NC_056532.1")
data$Chr = substr(data$Chr,1,11) #keep only the chromosome information
data = data[data$Chr %in% chrs,] #keep only the chromosomes of interest
data = data[complete.cases(data),]
rownames(data) = data[,1]
data = data[,-c(1:2)]
data = data[,-c(4,22)] #removes two bad merriami

spec.data = read.csv("morphological_data.csv")
spec.data = spec.data[spec.data$ID %in% colnames(data),]
group = paste(spec.data$Species, spec.data$Sex, spec.data$Treatment, sep = "")

DEGs = DGEList(counts = data, group = group)
design = model.matrix(~0+group, data = DEGs$samples)
keep = filterByExpr(DEGs, design = design, group = group)
DEGs = DEGs[keep,] #keeps 15234
data = data[rownames(data) %in% rownames(DEGs),]
DEGs = calcNormFactors(DEGs)
logcpm = cpm(DEGs, log = T)
DEGs = estimateDisp(DEGs, design)
DEGs$common.dispersion #0.181
fit_robust_auto = glmQLFit(DEGs, design, robust = T)

####undulatus vs virgatus####
#male
contrasts_auto = makeContrasts(und.TvC.male = groupSCUNMTEST - groupSCUNMCONT, levels = design)
QLF.und.TvC.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.TvC.male"])
de = decideTestsDGE(QLF.und.TvC.male, p = 0.05)
summary(de) #35 to test, 32 to controls

logfc.und.male = data.frame(topTags(QLF.und.TvC.male, n = Inf))
logfc.TEST.und.male = subset(logfc.und.male, logFC > 0 & FDR < 0.05)
logfc.CONT.und.male = subset(logfc.und.male, logFC < 0 & FDR < 0.05)

contrasts_auto = makeContrasts(virg.TvC.male = groupSCVIMTEST - groupSCVIMCONT, levels = design)
QLF.virg.TvC.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.TvC.male"])
de = decideTestsDGE(QLF.virg.TvC.male, p = 0.05)
summary(de) #343 to test, 90 to controls

logfc.virg.male = data.frame(topTags(QLF.virg.TvC.male, n = Inf))
logfc.TEST.virg.male = subset(logfc.virg.male, logFC > 0 & FDR < 0.05)
logfc.CONT.virg.male = subset(logfc.virg.male, logFC < 0 & FDR < 0.05)


logfc.und.male$gene = rownames(logfc.und.male)
logfc.virg.male$gene = rownames(logfc.virg.male)

data.merged = merge(logfc.und.male, logfc.virg.male, by = "gene")

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$gene))
colnames(fig.data) = c("undulatus", "virgatus", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$undulatus = as.numeric(fig.data$undulatus)
fig.data$virgatus = as.numeric(fig.data$virgatus)
fig.data$density = get_density(fig.data$undulatus, fig.data$virgatus, n = 1000)

#male main effect
contrasts_auto = makeContrasts(undandvirg.male = (groupSCUNMTEST + groupSCVIMTEST)/2 - (groupSCUNMCONT + groupSCVIMCONT)/2, levels = design)
QLF.undandvirg.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undandvirg.male"])
de = decideTestsDGE(QLF.undandvirg.male, p = 0.05)
summary(de) #168 down and 286 up
logfc.undandvirg.male = data.frame(topTags(QLF.undandvirg.male, n = Inf))
logfc.undandvirg.male.up = subset(logfc.undandvirg.male, logFC > 0 & FDR < 0.05)
logfc.undandvirg.male.down = subset(logfc.undandvirg.male, logFC < 0 & FDR < 0.05)

#male interaction
contrasts_auto = makeContrasts(undvsvirg.male = (groupSCUNMTEST - groupSCUNMCONT) - (groupSCVIMTEST - groupSCVIMCONT), levels = design)
QLF.undvsvirg.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undvsvirg.male"])
de = decideTestsDGE(QLF.undvsvirg.male, p = 0.05)
summary(de) #1 and 15
logfc.undvsvirg.male = data.frame(topTags(QLF.undvsvirg.male, n = Inf))
logfc.undvsvirg.male.up = subset(logfc.undvsvirg.male, logFC > 0 & FDR < 0.05)
logfc.undvsvirg.male.down = subset(logfc.undvsvirg.male, logFC < 0 & FDR < 0.05)

#count overlap
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandvirg.male.up[which(rownames(logfc.undandvirg.male.up) %in% rownames(logfc.undvsvirg.male.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandvirg.male.down[which(rownames(logfc.undandvirg.male.down) %in% rownames(logfc.undvsvirg.male.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandvirg.male.up[which(rownames(logfc.undandvirg.male.up) %in% rownames(logfc.undvsvirg.male.up)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandvirg.male.down[which(rownames(logfc.undandvirg.male.down) %in% rownames(logfc.undvsvirg.male.up)),]),])

und.virg.male = ggplot(fig.data, aes(x = undulatus, y = virgatus))+
  geom_point(aes(x = undulatus, y = virgatus, col = density))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.male.up),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.male.down),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsvirg.male.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsvirg.male.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.male.up[which(rownames(logfc.undandvirg.male.up) %in% rownames(logfc.undvsvirg.male.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.male.down[which(rownames(logfc.undandvirg.male.down) %in% rownames(logfc.undvsvirg.male.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.male.up[which(rownames(logfc.undandvirg.male.up) %in% rownames(logfc.undvsvirg.male.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.male.down[which(rownames(logfc.undandvirg.male.down) %in% rownames(logfc.undvsvirg.male.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  xlab(expression(bolditalic("S. undulatus ")*bold("Male log"["2"]*"FC")))+
  ylab(expression(bolditalic("S. virgatus ")*bold("Male log"["2"]*"FC")))+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6),expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$virgatus ~ fig.data$undulatus)) #r2 = 0.032, p < 0.001
cor.test(fig.data$undulatus, fig.data$virgatus)

#female
contrasts_auto = makeContrasts(und.TvC.female = groupSCUNFTEST - groupSCUNFCONT, levels = design)
QLF.und.TvC.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.TvC.female"])
de = decideTestsDGE(QLF.und.TvC.female, p = 0.05)
summary(de) #25 to test, 21 to controls

logfc.und.female = data.frame(topTags(QLF.und.TvC.female, n = Inf))
logfc.TEST.und.female = subset(logfc.und.female, logFC > 0 & FDR < 0.05)
logfc.CONT.und.female = subset(logfc.und.female, logFC < 0 & FDR < 0.05)

contrasts_auto = makeContrasts(virg.TvC.female = groupSCVIFTEST - groupSCVIFCONT, levels = design)
QLF.virg.TvC.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.TvC.female"])
de = decideTestsDGE(QLF.virg.TvC.female, p = 0.05)
summary(de) #99 to test, 21 to controls

logfc.virg.female = data.frame(topTags(QLF.virg.TvC.female, n = Inf))
logfc.TEST.virg.female = subset(logfc.virg.female, logFC > 0 & FDR < 0.05)
logfc.CONT.virg.female = subset(logfc.virg.female, logFC < 0 & FDR < 0.05)


logfc.und.female$gene = rownames(logfc.und.female)
logfc.virg.female$gene = rownames(logfc.virg.female)

data.merged = merge(logfc.und.female, logfc.virg.female, by = "gene")

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$gene))
colnames(fig.data) = c("undulatus", "virgatus", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$undulatus = as.numeric(fig.data$undulatus)
fig.data$virgatus = as.numeric(fig.data$virgatus)
fig.data$density = get_density(fig.data$undulatus, fig.data$virgatus, n = 1000)

#female main effect
contrasts_auto = makeContrasts(undandvirg.female = (groupSCUNFTEST + groupSCVIFTEST)/2 - (groupSCUNFCONT + groupSCVIFCONT)/2, levels = design)
QLF.undandvirg.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undandvirg.female"])
de = decideTestsDGE(QLF.undandvirg.female, p = 0.05)
summary(de) #59 down and 139 up
logfc.undandvirg.female = data.frame(topTags(QLF.undandvirg.female, n = Inf))
logfc.undandvirg.female.up = subset(logfc.undandvirg.female, logFC > 0 & FDR < 0.05)
logfc.undandvirg.female.down = subset(logfc.undandvirg.female, logFC < 0 & FDR < 0.05)

#female interaction
contrasts_auto = makeContrasts(undvsvirg.female = (groupSCUNFTEST - groupSCUNFCONT) - (groupSCVIFTEST - groupSCVIFCONT), levels = design)
QLF.undvsvirg.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undvsvirg.female"])
de = decideTestsDGE(QLF.undvsvirg.female, p = 0.05)
summary(de) #1 and 0
logfc.undvsvirg.female = data.frame(topTags(QLF.undvsvirg.female, n = Inf))
logfc.undvsvirg.female.up = subset(logfc.undvsvirg.female, logFC > 0 & FDR < 0.05)
logfc.undvsvirg.female.down = subset(logfc.undvsvirg.female, logFC < 0 & FDR < 0.05)

und.virg.female = ggplot(fig.data, aes(x = undulatus, y = virgatus))+
  geom_point(aes(x = undulatus, y = virgatus), col = "grey")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.female.up),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.female.down),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsvirg.female.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsvirg.female.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.female.up[which(rownames(logfc.undandvirg.female.up) %in% rownames(logfc.undvsvirg.female.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+ 
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.female.down[which(rownames(logfc.undandvirg.female.down) %in% rownames(logfc.undvsvirg.female.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.female.up[which(rownames(logfc.undandvirg.female.up) %in% rownames(logfc.undvsvirg.female.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.female.down[which(rownames(logfc.undandvirg.female.down) %in% rownames(logfc.undvsvirg.female.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  xlab(expression(bolditalic("S. undulatus ")*bold("Female log"["2"]*"FC")))+
  ylab(expression(bolditalic("S. virgatus ")*bold("Female log"["2"]*"FC")))+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$virgatus ~ fig.data$undulatus)) #r2 = 0.017

cor.test(fig.data$undulatus, fig.data$virgatus)

#both
contrasts_auto = makeContrasts(und.TvC.both = (groupSCUNFTEST + groupSCUNMTEST)/2 - (groupSCUNFCONT + groupSCUNMCONT)/2, levels = design)
QLF.und.TvC.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.TvC.both"])
de = decideTestsDGE(QLF.und.TvC.both, p = 0.05)
summary(de) #146 to test, 167 to controls

logfc.und.both = data.frame(topTags(QLF.und.TvC.both, n = Inf))
logfc.TEST.und.both = subset(logfc.und.both, logFC > 0 & FDR < 0.05)
logfc.CONT.und.both = subset(logfc.und.both, logFC < 0 & FDR < 0.05)

contrasts_auto = makeContrasts(virg.TvC.both = (groupSCVIFTEST + groupSCVIMTEST)/2 - (groupSCVIFCONT + groupSCVIMCONT)/2, levels = design)
QLF.virg.TvC.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.TvC.both"])
de = decideTestsDGE(QLF.virg.TvC.both, p = 0.05)
summary(de) #704 to test, 222 to controls

logfc.virg.both = data.frame(topTags(QLF.virg.TvC.both, n = Inf))
logfc.TEST.virg.both = subset(logfc.virg.both, logFC > 0 & FDR < 0.05)
logfc.CONT.virg.both = subset(logfc.virg.both, logFC < 0 & FDR < 0.05)

logfc.und.both$gene = rownames(logfc.und.both)
logfc.virg.both$gene = rownames(logfc.virg.both)

data.merged = merge(logfc.und.both, logfc.virg.both, by = "gene")

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$gene))
colnames(fig.data) = c("undulatus", "virgatus", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$undulatus = as.numeric(fig.data$undulatus)
fig.data$virgatus = as.numeric(fig.data$virgatus)
fig.data$density = get_density(fig.data$undulatus, fig.data$virgatus, n = 1000)

#both main effect
contrasts_auto = makeContrasts(undandvirg.both = (groupSCUNFTEST + groupSCVIFTEST + groupSCUNMTEST + groupSCVIMTEST)/4 - (groupSCUNFCONT + groupSCVIFCONT + groupSCUNMCONT + groupSCVIMCONT)/4, levels = design)
QLF.undandvirg.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undandvirg.both"])
de = decideTestsDGE(QLF.undandvirg.both, p = 0.05)
summary(de) #539 down and 743 up
logfc.undandvirg.both = data.frame(topTags(QLF.undandvirg.both, n = Inf))
logfc.undandvirg.both.up = subset(logfc.undandvirg.both, logFC > 0 & FDR < 0.05)
logfc.undandvirg.both.down = subset(logfc.undandvirg.both, logFC < 0 & FDR < 0.05)

#both interaction
contrasts_auto = makeContrasts(undvsvirg.both = ((groupSCUNFTEST + groupSCUNMTEST)/2 - (groupSCUNFCONT + groupSCUNMCONT)/2) - ((groupSCVIFTEST + groupSCVIMTEST)/2 - (groupSCVIFCONT + groupSCVIMCONT)/2), levels = design)
QLF.undvsvirg.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undvsvirg.both"])
de = decideTestsDGE(QLF.undvsvirg.both, p = 0.05)
summary(de) #83 and 9
logfc.undvsvirg.both = data.frame(topTags(QLF.undvsvirg.both, n = Inf))
logfc.undvsvirg.both.up = subset(logfc.undvsvirg.both, logFC > 0 & FDR < 0.05)
logfc.undvsvirg.both.down = subset(logfc.undvsvirg.both, logFC < 0 & FDR < 0.05)

nrow(data.merged[data.merged$gene %in% rownames(logfc.undandvirg.both.up[which(rownames(logfc.undandvirg.both.up) %in% rownames(logfc.undvsvirg.both.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandvirg.both.down[which(rownames(logfc.undandvirg.both.down) %in% rownames(logfc.undvsvirg.both.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandvirg.both.up[which(rownames(logfc.undandvirg.both.up) %in% rownames(logfc.undvsvirg.both.up)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandvirg.both.down[which(rownames(logfc.undandvirg.both.down) %in% rownames(logfc.undvsvirg.both.up)),]),])

und.virg.both = ggplot(fig.data, aes(x = undulatus, y = virgatus))+
  geom_point(aes(x = undulatus, y = virgatus), col = "grey")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.both.up),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.both.down),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsvirg.both.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsvirg.both.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.both.up[which(rownames(logfc.undandvirg.both.up) %in% rownames(logfc.undvsvirg.both.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.both.down[which(rownames(logfc.undandvirg.both.down) %in% rownames(logfc.undvsvirg.both.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.both.up[which(rownames(logfc.undandvirg.both.up) %in% rownames(logfc.undvsvirg.both.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandvirg.both.down[which(rownames(logfc.undandvirg.both.down) %in% rownames(logfc.undvsvirg.both.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  xlab(expression(bolditalic("S. undulatus ")*bold("log"["2"]*"FC")))+
  ylab(expression(bolditalic("S. virgatus ")*bold("log"["2"]*"FC")))+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$virgatus ~ fig.data$undulatus)) #r2 = 0.081
cor.test(fig.data$undulatus, fig.data$virgatus)


####undulatus vs merriami####
#male
contrasts_auto = makeContrasts(und.TvC.male = groupSCUNMTEST - groupSCUNMCONT, levels = design)
QLF.und.TvC.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.TvC.male"])
de = decideTestsDGE(QLF.und.TvC.male, p = 0.05)
summary(de) #35 to test, 32 to controls

logfc.und.male = data.frame(topTags(QLF.und.TvC.male, n = Inf))
logfc.TEST.und.male = subset(logfc.und.male, logFC > 0 & FDR < 0.05)
logfc.CONT.und.male = subset(logfc.und.male, logFC < 0 & FDR < 0.05)

contrasts_auto = makeContrasts(mer.TvC.male = groupSCMEMTEST - groupSCMEMCONT, levels = design)
QLF.mer.TvC.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.TvC.male"])
de = decideTestsDGE(QLF.mer.TvC.male, p = 0.05)
summary(de) #316 to test, 229 to controls

logfc.mer.male = data.frame(topTags(QLF.mer.TvC.male, n = Inf))
logfc.TEST.mer.male = subset(logfc.mer.male, logFC > 0 & FDR < 0.05)
logfc.CONT.mer.male = subset(logfc.mer.male, logFC < 0 & FDR < 0.05)


logfc.und.male$gene = rownames(logfc.und.male)
logfc.mer.male$gene = rownames(logfc.mer.male)

data.merged = merge(logfc.und.male, logfc.mer.male, by = "gene")

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$gene))
colnames(fig.data) = c("undulatus", "merriami", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$undulatus = as.numeric(fig.data$undulatus)
fig.data$merriami = as.numeric(fig.data$merriami)
fig.data$density = get_density(fig.data$undulatus, fig.data$merriami, n = 1000)

#male main effect
contrasts_auto = makeContrasts(undandmer.male = (groupSCUNMTEST + groupSCMEMTEST)/2 - (groupSCUNMCONT + groupSCMEMCONT)/2, levels = design)
QLF.undandmer.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undandmer.male"])
de = decideTestsDGE(QLF.undandmer.male, p = 0.05)
summary(de) #321 down and 296 up
logfc.undandmer.male = data.frame(topTags(QLF.undandmer.male, n = Inf))
logfc.undandmer.male.up = subset(logfc.undandmer.male, logFC > 0 & FDR < 0.05)
logfc.undandmer.male.down = subset(logfc.undandmer.male, logFC < 0 & FDR < 0.05)

#male interaction
contrasts_auto = makeContrasts(undvsmer.male = (groupSCUNMTEST - groupSCUNMCONT) - (groupSCMEMTEST - groupSCMEMCONT), levels = design)
QLF.undvsmer.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undvsmer.male"])
de = decideTestsDGE(QLF.undvsmer.male, p = 0.05)
summary(de) #27 and 61
logfc.undvsmer.male = data.frame(topTags(QLF.undvsmer.male, n = Inf))
logfc.undvsmer.male.up = subset(logfc.undvsmer.male, logFC > 0 & FDR < 0.05)
logfc.undvsmer.male.down = subset(logfc.undvsmer.male, logFC < 0 & FDR < 0.05)

nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.male.up[which(rownames(logfc.undandmer.male.up) %in% rownames(logfc.undvsmer.male.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.male.down[which(rownames(logfc.undandmer.male.down) %in% rownames(logfc.undvsmer.male.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.male.up[which(rownames(logfc.undandmer.male.up) %in% rownames(logfc.undvsmer.male.up)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.male.down[which(rownames(logfc.undandmer.male.down) %in% rownames(logfc.undvsmer.male.up)),]),])

und.mer.male = ggplot(fig.data, aes(x = undulatus, y = merriami))+
  geom_point(aes(x = undulatus, y = merriami, col = density))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.male.up),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.male.down),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsmer.male.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsmer.male.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.male.up[which(rownames(logfc.undandmer.male.up) %in% rownames(logfc.undvsmer.male.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+ 
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.male.down[which(rownames(logfc.undandmer.male.down) %in% rownames(logfc.undvsmer.male.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.male.up[which(rownames(logfc.undandmer.male.up) %in% rownames(logfc.undvsmer.male.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.male.down[which(rownames(logfc.undandmer.male.down) %in% rownames(logfc.undvsmer.male.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  xlab(expression(bolditalic("S. undulatus ")*bold("Male log"["2"]*"FC")))+
  ylab(expression(bolditalic("S. merrami ")*bold("Male log"["2"]*"FC")))+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$merriami ~ fig.data$undulatus)) #r2 = 0.030
cor.test(fig.data$undulatus, fig.data$merriami)

#female
contrasts_auto = makeContrasts(und.TvC.female = groupSCUNFTEST - groupSCUNFCONT, levels = design)
QLF.und.TvC.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.TvC.female"])
de = decideTestsDGE(QLF.und.TvC.female, p = 0.05)
summary(de) #25 to test, 21 to controls

logfc.und.female = data.frame(topTags(QLF.und.TvC.female, n = Inf))
logfc.TEST.und.female = subset(logfc.und.female, logFC > 0 & FDR < 0.05)
logfc.CONT.und.female = subset(logfc.und.female, logFC < 0 & FDR < 0.05)

contrasts_auto = makeContrasts(mer.TvC.female = groupSCMEFTEST - groupSCMEFCONT, levels = design)
QLF.mer.TvC.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.TvC.female"])
de = decideTestsDGE(QLF.mer.TvC.female, p = 0.05)
summary(de) #168 to test, 98 to controls

logfc.mer.female = data.frame(topTags(QLF.mer.TvC.female, n = Inf))
logfc.TEST.mer.female = subset(logfc.mer.female, logFC > 0 & FDR < 0.05)
logfc.CONT.mer.female = subset(logfc.mer.female, logFC < 0 & FDR < 0.05)


logfc.und.female$gene = rownames(logfc.und.female)
logfc.mer.female$gene = rownames(logfc.mer.female)

data.merged = merge(logfc.und.female, logfc.mer.female, by = "gene")

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$gene))
colnames(fig.data) = c("undulatus", "merriami", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$undulatus = as.numeric(fig.data$undulatus)
fig.data$merriami = as.numeric(fig.data$merriami)
fig.data$density = get_density(fig.data$undulatus, fig.data$merriami, n = 1000)

#female main effect
contrasts_auto = makeContrasts(undandmer.female = (groupSCUNFTEST + groupSCMEFTEST)/2 - (groupSCUNFCONT + groupSCMEFCONT)/2, levels = design)
QLF.undandmer.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undandmer.female"])
de = decideTestsDGE(QLF.undandmer.female, p = 0.05)
summary(de) #51 down and 82 up
logfc.undandmer.female = data.frame(topTags(QLF.undandmer.female, n = Inf))
logfc.undandmer.female.up = subset(logfc.undandmer.female, logFC > 0 & FDR < 0.05)
logfc.undandmer.female.down = subset(logfc.undandmer.female, logFC < 0 & FDR < 0.05)

#female interaction
contrasts_auto = makeContrasts(undvsmer.female = (groupSCUNFTEST - groupSCUNFCONT) - (groupSCMEFTEST - groupSCMEFCONT), levels = design)
QLF.undvsmer.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undvsmer.female"])
de = decideTestsDGE(QLF.undvsmer.female, p = 0.05)
summary(de) #6 and 23
logfc.undvsmer.female = data.frame(topTags(QLF.undvsmer.female, n = Inf))
logfc.undvsmer.female.up = subset(logfc.undvsmer.female, logFC > 0 & FDR < 0.05)
logfc.undvsmer.female.down = subset(logfc.undvsmer.female, logFC < 0 & FDR < 0.05)

nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.female.up[which(rownames(logfc.undandmer.female.up) %in% rownames(logfc.undvsmer.female.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.female.down[which(rownames(logfc.undandmer.female.down) %in% rownames(logfc.undvsmer.female.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.female.up[which(rownames(logfc.undandmer.female.up) %in% rownames(logfc.undvsmer.female.up)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.female.down[which(rownames(logfc.undandmer.female.down) %in% rownames(logfc.undvsmer.female.up)),]),])

und.mer.female = ggplot(fig.data, aes(x = undulatus, y = merriami))+
  geom_point(aes(x = undulatus, y = merriami, col = density))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.female.up),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.female.down),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsmer.female.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsmer.female.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.female.up[which(rownames(logfc.undandmer.female.up) %in% rownames(logfc.undvsmer.female.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+ 
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.female.down[which(rownames(logfc.undandmer.female.down) %in% rownames(logfc.undvsmer.female.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.female.up[which(rownames(logfc.undandmer.female.up) %in% rownames(logfc.undvsmer.female.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.female.down[which(rownames(logfc.undandmer.female.down) %in% rownames(logfc.undvsmer.female.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  xlab(expression(bolditalic("S. undulatus ")*bold("Female log"["2"]*"FC")))+
  ylab(expression(bolditalic("S. merrami ")*bold("Female log"["2"]*"FC")))+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$merriami ~ fig.data$undulatus)) #r2 = 0.020
cor.test(fig.data$undulatus, fig.data$merriami)

#both
contrasts_auto = makeContrasts(und.TvC.both = (groupSCUNFTEST + groupSCUNMTEST)/2 - (groupSCUNFCONT + groupSCUNMCONT)/2, levels = design)
QLF.und.TvC.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"und.TvC.both"])
de = decideTestsDGE(QLF.und.TvC.both, p = 0.05)
summary(de) #146 to test, 167 to controls

logfc.und.both = data.frame(topTags(QLF.und.TvC.both, n = Inf))
logfc.TEST.und.both = subset(logfc.und.both, logFC > 0 & FDR < 0.05)
logfc.CONT.und.both = subset(logfc.und.both, logFC < 0 & FDR < 0.05)

contrasts_auto = makeContrasts(mer.TvC.both = (groupSCMEFTEST + groupSCMEMTEST)/2 - (groupSCMEFCONT + groupSCMEMCONT)/2, levels = design)
QLF.mer.TvC.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.TvC.both"])
de = decideTestsDGE(QLF.mer.TvC.both, p = 0.05)
summary(de) #573 to test, 444 to controls

logfc.mer.both = data.frame(topTags(QLF.mer.TvC.both, n = Inf))
logfc.TEST.mer.both = subset(logfc.mer.both, logFC > 0 & FDR < 0.05)
logfc.CONT.mer.both = subset(logfc.mer.both, logFC < 0 & FDR < 0.05)


logfc.und.both$gene = rownames(logfc.und.both)
logfc.mer.both$gene = rownames(logfc.mer.both)

data.merged = merge(logfc.und.both, logfc.mer.both, by = "gene")

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$gene))
colnames(fig.data) = c("undulatus", "merriami", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$undulatus = as.numeric(fig.data$undulatus)
fig.data$merriami = as.numeric(fig.data$merriami)
fig.data$density = get_density(fig.data$undulatus, fig.data$merriami, n = 1000)

#both main effect
contrasts_auto = makeContrasts(undandmer.both = (groupSCUNFTEST + groupSCMEFTEST + groupSCUNMTEST + groupSCMEMTEST)/4 - (groupSCUNFCONT + groupSCMEFCONT + groupSCUNMCONT + groupSCMEMCONT)/4, levels = design)
QLF.undandmer.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undandmer.both"])
de = decideTestsDGE(QLF.undandmer.both, p = 0.05)
summary(de) #549 down and 557 up
logfc.undandmer.both = data.frame(topTags(QLF.undandmer.both, n = Inf))
logfc.undandmer.both.up = subset(logfc.undandmer.both, logFC > 0 & FDR < 0.05)
logfc.undandmer.both.down = subset(logfc.undandmer.both, logFC < 0 & FDR < 0.05)

#both interaction
contrasts_auto = makeContrasts(undvsmer.both = ((groupSCUNFTEST + groupSCUNMTEST)/2 - (groupSCUNFCONT + groupSCUNMCONT)/2) - ((groupSCMEFTEST + groupSCMEMTEST)/2 - (groupSCMEFCONT + groupSCMEMCONT)/2), levels = design)
QLF.undvsmer.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"undvsmer.both"])
de = decideTestsDGE(QLF.undvsmer.both, p = 0.05)
summary(de) #141 and 234
logfc.undvsmer.both = data.frame(topTags(QLF.undvsmer.both, n = Inf))
logfc.undvsmer.both.up = subset(logfc.undvsmer.both, logFC > 0 & FDR < 0.05)
logfc.undvsmer.both.down = subset(logfc.undvsmer.both, logFC < 0 & FDR < 0.05)

nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.both.up[which(rownames(logfc.undandmer.both.up) %in% rownames(logfc.undvsmer.both.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.both.down[which(rownames(logfc.undandmer.both.down) %in% rownames(logfc.undvsmer.both.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.both.up[which(rownames(logfc.undandmer.both.up) %in% rownames(logfc.undvsmer.both.up)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.undandmer.both.down[which(rownames(logfc.undandmer.both.down) %in% rownames(logfc.undvsmer.both.up)),]),])

und.mer.both = ggplot(fig.data, aes(x = undulatus, y = merriami))+
  geom_point(aes(x = undulatus, y = merriami, col = density))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.both.up),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.both.down),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsmer.both.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undvsmer.both.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.both.up[which(rownames(logfc.undandmer.both.up) %in% rownames(logfc.undvsmer.both.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+ 
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.both.down[which(rownames(logfc.undandmer.both.down) %in% rownames(logfc.undvsmer.both.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.both.up[which(rownames(logfc.undandmer.both.up) %in% rownames(logfc.undvsmer.both.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.undandmer.both.down[which(rownames(logfc.undandmer.both.down) %in% rownames(logfc.undvsmer.both.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  xlab(expression(bolditalic("S. undulatus ")*bold("log"["2"]*"FC")))+
  ylab(expression(bolditalic("S. merrami ")*bold("log"["2"]*"FC")))+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$merriami ~ fig.data$undulatus)) #r2 = 0.041
cor.test(fig.data$undulatus, fig.data$merriami)


####virgatus vs merriami####
#male
contrasts_auto = makeContrasts(virg.TvC.male = groupSCVIMTEST - groupSCVIMCONT, levels = design)
QLF.virg.TvC.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.TvC.male"])
de = decideTestsDGE(QLF.virg.TvC.male, p = 0.05)
summary(de) #343 to test, 90 to controls

logfc.virg.male = data.frame(topTags(QLF.virg.TvC.male, n = Inf))
logfc.TEST.virg.male = subset(logfc.virg.male, logFC > 0 & FDR < 0.05)
logfc.CONT.virg.male = subset(logfc.virg.male, logFC < 0 & FDR < 0.05)

contrasts_auto = makeContrasts(mer.TvC.male = groupSCMEMTEST - groupSCMEMCONT, levels = design)
QLF.mer.TvC.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.TvC.male"])
de = decideTestsDGE(QLF.mer.TvC.male, p = 0.05)
summary(de) #316 to test, 229 to controls

logfc.mer.male = data.frame(topTags(QLF.mer.TvC.male, n = Inf))
logfc.TEST.mer.male = subset(logfc.mer.male, logFC > 0 & FDR < 0.05)
logfc.CONT.mer.male = subset(logfc.mer.male, logFC < 0 & FDR < 0.05)


logfc.virg.male$gene = rownames(logfc.virg.male)
logfc.mer.male$gene = rownames(logfc.mer.male)

data.merged = merge(logfc.virg.male, logfc.mer.male, by = "gene")

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$gene))
colnames(fig.data) = c("virgatus", "merriami", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$virgatus = as.numeric(fig.data$virgatus)
fig.data$merriami = as.numeric(fig.data$merriami)
fig.data$density = get_density(fig.data$virgatus, fig.data$merriami, n = 1000)

#male main effect
contrasts_auto = makeContrasts(virgandmer.male = (groupSCVIMTEST + groupSCMEMTEST)/2 - (groupSCVIMCONT + groupSCMEMCONT)/2, levels = design)
QLF.virgandmer.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virgandmer.male"])
de = decideTestsDGE(QLF.virgandmer.male, p = 0.05)
summary(de) #271 down and 502 up
logfc.virgandmer.male = data.frame(topTags(QLF.virgandmer.male, n = Inf))
logfc.virgandmer.male.up = subset(logfc.virgandmer.male, logFC > 0 & FDR < 0.05)
logfc.virgandmer.male.down = subset(logfc.virgandmer.male, logFC < 0 & FDR < 0.05)

#male interaction
contrasts_auto = makeContrasts(virgvsmer.male = (groupSCVIMTEST - groupSCVIMCONT) - (groupSCMEMTEST - groupSCMEMCONT), levels = design)
QLF.virgvsmer.male = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virgvsmer.male"])
de = decideTestsDGE(QLF.virgvsmer.male, p = 0.05)
summary(de) #153 and 102
logfc.virgvsmer.male = data.frame(topTags(QLF.virgvsmer.male, n = Inf))
logfc.virgvsmer.male.up = subset(logfc.virgvsmer.male, logFC > 0 & FDR < 0.05)
logfc.virgvsmer.male.down = subset(logfc.virgvsmer.male, logFC < 0 & FDR < 0.05)

nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.male.up[which(rownames(logfc.virgandmer.male.up) %in% rownames(logfc.virgvsmer.male.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.male.down[which(rownames(logfc.virgandmer.male.down) %in% rownames(logfc.virgvsmer.male.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.male.up[which(rownames(logfc.virgandmer.male.up) %in% rownames(logfc.virgvsmer.male.up)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.male.down[which(rownames(logfc.virgandmer.male.down) %in% rownames(logfc.virgvsmer.male.up)),]),])

virg.mer.male = ggplot(fig.data, aes(x = virgatus, y = merriami))+
  geom_point(aes(x = virgatus, y = merriami, col = density))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.male.up),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.male.down),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgvsmer.male.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgvsmer.male.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.male.up[which(rownames(logfc.virgandmer.male.up) %in% rownames(logfc.virgvsmer.male.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+ 
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.male.down[which(rownames(logfc.virgandmer.male.down) %in% rownames(logfc.virgvsmer.male.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.male.up[which(rownames(logfc.virgandmer.male.up) %in% rownames(logfc.virgvsmer.male.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.male.down[which(rownames(logfc.virgandmer.male.down) %in% rownames(logfc.virgvsmer.male.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  xlab(expression(bolditalic("S. virgatus ")*bold("Male log"["2"]*"FC")))+
  ylab(expression(bolditalic("S. merrami ")*bold("Male log"["2"]*"FC")))+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$merriami ~ fig.data$virgatus)) #r2 = 0.024
cor.test(fig.data$virgatus, fig.data$merriami)

#female
contrasts_auto = makeContrasts(virg.TvC.female = groupSCVIFTEST - groupSCVIFCONT, levels = design)
QLF.virg.TvC.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.TvC.female"])
de = decideTestsDGE(QLF.virg.TvC.female, p = 0.05)
summary(de) #99 to test, 21 to controls

logfc.virg.female = data.frame(topTags(QLF.virg.TvC.female, n = Inf))
logfc.TEST.virg.female = subset(logfc.virg.female, logFC > 0 & FDR < 0.05)
logfc.CONT.virg.female = subset(logfc.virg.female, logFC < 0 & FDR < 0.05)

contrasts_auto = makeContrasts(mer.TvC.female = groupSCMEFTEST - groupSCMEFCONT, levels = design)
QLF.mer.TvC.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.TvC.female"])
de = decideTestsDGE(QLF.mer.TvC.female, p = 0.05)
summary(de) #168 to test, 98 to controls

logfc.mer.female = data.frame(topTags(QLF.mer.TvC.female, n = Inf))
logfc.TEST.mer.female = subset(logfc.mer.female, logFC > 0 & FDR < 0.05)
logfc.CONT.mer.female = subset(logfc.mer.female, logFC < 0 & FDR < 0.05)


logfc.virg.female$gene = rownames(logfc.virg.female)
logfc.mer.female$gene = rownames(logfc.mer.female)

data.merged = merge(logfc.virg.female, logfc.mer.female, by = "gene")

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$gene))
colnames(fig.data) = c("virgatus", "merriami", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$virgatus = as.numeric(fig.data$virgatus)
fig.data$merriami = as.numeric(fig.data$merriami)
fig.data$density = get_density(fig.data$virgatus, fig.data$merriami, n = 1000)

#female main effect
contrasts_auto = makeContrasts(virgandmer.female = (groupSCVIFTEST + groupSCMEFTEST)/2 - (groupSCVIFCONT + groupSCMEFCONT)/2, levels = design)
QLF.virgandmer.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virgandmer.female"])
de = decideTestsDGE(QLF.virgandmer.female, p = 0.05)
summary(de) #48 down and 139 up
logfc.virgandmer.female = data.frame(topTags(QLF.virgandmer.female, n = Inf))
logfc.virgandmer.female.up = subset(logfc.virgandmer.female, logFC > 0 & FDR < 0.05)
logfc.virgandmer.female.down = subset(logfc.virgandmer.female, logFC < 0 & FDR < 0.05)

#female interaction
contrasts_auto = makeContrasts(virgvsmer.female = (groupSCVIFTEST - groupSCVIFCONT) - (groupSCMEFTEST - groupSCMEFCONT), levels = design)
QLF.virgvsmer.female = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virgvsmer.female"])
de = decideTestsDGE(QLF.virgvsmer.female, p = 0.05)
summary(de) #12 and 13
logfc.virgvsmer.female = data.frame(topTags(QLF.virgvsmer.female, n = Inf))
logfc.virgvsmer.female.up = subset(logfc.virgvsmer.female, logFC > 0 & FDR < 0.05)
logfc.virgvsmer.female.down = subset(logfc.virgvsmer.female, logFC < 0 & FDR < 0.05)

nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.female.up[which(rownames(logfc.virgandmer.female.up) %in% rownames(logfc.virgvsmer.female.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.female.down[which(rownames(logfc.virgandmer.female.down) %in% rownames(logfc.virgvsmer.female.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.female.up[which(rownames(logfc.virgandmer.female.up) %in% rownames(logfc.virgvsmer.female.up)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.female.down[which(rownames(logfc.virgandmer.female.down) %in% rownames(logfc.virgvsmer.female.up)),]),])

virg.mer.female = ggplot(fig.data, aes(x = virgatus, y = merriami))+
  geom_point(aes(x = virgatus, y = merriami, col = density))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.female.up),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.female.down),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgvsmer.female.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgvsmer.female.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.female.up[which(rownames(logfc.virgandmer.female.up) %in% rownames(logfc.virgvsmer.female.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+ 
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.female.down[which(rownames(logfc.virgandmer.female.down) %in% rownames(logfc.virgvsmer.female.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.female.up[which(rownames(logfc.virgandmer.female.up) %in% rownames(logfc.virgvsmer.female.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.female.down[which(rownames(logfc.virgandmer.female.down) %in% rownames(logfc.virgvsmer.female.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  xlab(expression(bolditalic("S. undulatus ")*bold("Female log"["2"]*"FC")))+
  ylab(expression(bolditalic("S. merrami ")*bold("Female log"["2"]*"FC")))+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$merriami ~ fig.data$virgatus)) #r2 = 0.0002
cor.test(fig.data$virgatus, fig.data$merriami)

#both
contrasts_auto = makeContrasts(virg.TvC.both = (groupSCVIFTEST + groupSCVIMTEST)/2 - (groupSCVIFCONT + groupSCVIMCONT)/2, levels = design)
QLF.virg.TvC.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virg.TvC.both"])
de = decideTestsDGE(QLF.virg.TvC.both, p = 0.05)
summary(de) #704 to test, 222 to controls

logfc.virg.both = data.frame(topTags(QLF.virg.TvC.both, n = Inf))
logfc.TEST.virg.both = subset(logfc.virg.both, logFC > 0 & FDR < 0.05)
logfc.CONT.virg.both = subset(logfc.virg.both, logFC < 0 & FDR < 0.05)

contrasts_auto = makeContrasts(mer.TvC.both = (groupSCMEFTEST + groupSCMEMTEST)/2 - (groupSCMEFCONT + groupSCMEMCONT)/2, levels = design)
QLF.mer.TvC.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"mer.TvC.both"])
de = decideTestsDGE(QLF.mer.TvC.both, p = 0.05)
summary(de) #573 to test, 444 to controls

logfc.mer.both = data.frame(topTags(QLF.mer.TvC.both, n = Inf))
logfc.TEST.mer.both = subset(logfc.mer.both, logFC > 0 & FDR < 0.05)
logfc.CONT.mer.both = subset(logfc.mer.both, logFC < 0 & FDR < 0.05)


logfc.virg.both$gene = rownames(logfc.virg.both)
logfc.mer.both$gene = rownames(logfc.mer.both)

data.merged = merge(logfc.virg.both, logfc.mer.both, by = "gene")

fig.data = data.frame(cbind(data.merged$logFC.x, data.merged$logFC.y, data.merged$gene))
colnames(fig.data) = c("virgatus", "merriami", "gene")
rownames(fig.data) = fig.data$gene
fig.data = fig.data[,-3]
fig.data$virgatus = as.numeric(fig.data$virgatus)
fig.data$merriami = as.numeric(fig.data$merriami)
fig.data$density = get_density(fig.data$virgatus, fig.data$merriami, n = 1000)

#both main effect
contrasts_auto = makeContrasts(virgandmer.both = (groupSCVIFTEST + groupSCMEFTEST + groupSCVIMTEST + groupSCMEMTEST)/4 - (groupSCVIFCONT + groupSCMEFCONT + groupSCVIMCONT + groupSCMEMCONT)/4, levels = design)
QLF.virgandmer.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virgandmer.both"])
de = decideTestsDGE(QLF.virgandmer.both, p = 0.05)
summary(de) #517 down and 875 up
logfc.virgandmer.both = data.frame(topTags(QLF.virgandmer.both, n = Inf))
logfc.virgandmer.both.up = subset(logfc.virgandmer.both, logFC > 0 & FDR < 0.05)
logfc.virgandmer.both.down = subset(logfc.virgandmer.both, logFC < 0 & FDR < 0.05)

#both interaction
contrasts_auto = makeContrasts(virgvsmer.both = ((groupSCVIFTEST + groupSCVIMTEST)/2 - (groupSCVIFCONT + groupSCVIMCONT)/2) - ((groupSCMEFTEST + groupSCMEMTEST)/2 - (groupSCMEFCONT + groupSCMEMCONT)/2), levels = design)
QLF.virgvsmer.both = glmQLFTest(fit_robust_auto, contrast = contrasts_auto[,"virgvsmer.both"])
de = decideTestsDGE(QLF.virgvsmer.both, p = 0.05)
summary(de) #333 and 246
logfc.virgvsmer.both = data.frame(topTags(QLF.virgvsmer.both, n = Inf))
logfc.virgvsmer.both.up = subset(logfc.virgvsmer.both, logFC > 0 & FDR < 0.05)
logfc.virgvsmer.both.down = subset(logfc.virgvsmer.both, logFC < 0 & FDR < 0.05)

nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.both.up[which(rownames(logfc.virgandmer.both.up) %in% rownames(logfc.virgvsmer.both.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.both.down[which(rownames(logfc.virgandmer.both.down) %in% rownames(logfc.virgvsmer.both.down)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.both.up[which(rownames(logfc.virgandmer.both.up) %in% rownames(logfc.virgvsmer.both.up)),]),])
nrow(data.merged[data.merged$gene %in% rownames(logfc.virgandmer.both.down[which(rownames(logfc.virgandmer.both.down) %in% rownames(logfc.virgvsmer.both.up)),]),])

virg.mer.both = ggplot(fig.data, aes(x = virgatus, y = merriami))+
  geom_point(aes(x = virgatus, y = merriami, col = density))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.both.up),],
             aes(x = logFC.x, y = logFC.y), fill = "dodgerblue3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.both.down),],
             aes(x = logFC.x, y = logFC.y), fill = "indianred1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgvsmer.both.up),],
             aes(x = logFC.x, y = logFC.y), fill = "darkgreen", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgvsmer.both.down),],
             aes(x = logFC.x, y = logFC.y), fill = "goldenrod1", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.both.up[which(rownames(logfc.virgandmer.both.up) %in% rownames(logfc.virgvsmer.both.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+ 
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.both.down[which(rownames(logfc.virgandmer.both.down) %in% rownames(logfc.virgvsmer.both.down)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.both.up[which(rownames(logfc.virgandmer.both.up) %in% rownames(logfc.virgvsmer.both.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  geom_point(data = data.merged[data.merged$gene %in% rownames(logfc.virgandmer.both.down[which(rownames(logfc.virgandmer.both.down) %in% rownames(logfc.virgvsmer.both.up)),]),],
             aes(x = logFC.x, y = logFC.y), fill = "green3", size = 2.5, pch = 21)+
  xlab(expression(bolditalic("S. undulatus ")*bold("log"["2"]*"FC")))+
  ylab(expression(bolditalic("S. merrami ")*bold("log"["2"]*"FC")))+
  scale_color_gradient(low = "grey75", high = "grey25")+
  coord_cartesian(xlim = c(-6,6), ylim = c(-6,6), expand = F)+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6))+
  geom_smooth(method = "lm", col = "black", se = F, lwd = 2)+
  plot_theme
summary(lm(fig.data$merriami ~ fig.data$virgatus)) #r2 = 0.020
cor.test(fig.data$virgatus, fig.data$merriami)


####testing differences in proportions####
#undvvirg v undvmer
prop.test(c(92,375), c(1328, 1310))
#x2 = 211.63, p < 0.001

#undvvirg v virgvmer
prop.test(c(92, 579), c(1328, 1752))
#x2 = 300.94, p < 0.001

#undvmer vs virgvmer
prop.test(c(375, 579), c(1310, 1752))
#x2 = 6.63, p = 0.01

####Testosterone figure####
spec.data = read.csv("morphological_data.csv")
spec.data = spec.data[spec.data$ID %in% colnames(data),]

spec.data$Species = factor(spec.data$Species, levels = c("SCUN", "SCVI", "SCME"))

test.fig = ggplot(spec.data, aes(x = Species, y = Test_conc, fill = Treatment))+
  geom_boxplot(width = 0.5, lwd = 1, col = "black", position = position_dodge(0.75),
               outlier.colour = "white")+
  scale_fill_manual(values=c("indianred1","indianred1","dodgerblue3","dodgerblue3"))+
  geom_point(size = 5, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25),
             pch = ifelse(spec.data$Sex == "M", 21, 24),
             stroke = 1.5,
             aes(x = Species, y = Test_conc, fill = ifelse(Treatment == "TEST", "indianred1","dodgerblue3")))+
  scale_y_continuous(breaks = seq(0,75,by=25))+
  coord_cartesian(ylim = c(-2,75), expand = F)+
  scale_x_discrete(labels = c("S. undulatus", "S. virgatus", "S. merriami"))+
  ylab("Plasma [Testosterone] (ng/ml)")+xlab("Species")+
  plot_theme

summary(aov(Test_conc ~ Species * Sex * Treatment, data = spec.data)) #only a treatment effect