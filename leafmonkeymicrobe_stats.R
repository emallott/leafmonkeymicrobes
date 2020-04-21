#Alpha diversity----
library(tidyverse)
library(nlme)
library(car)
library(multcomp)
library(lme4)
library(performance)

metadata = read.csv("leaf_map_r_8000.csv", header = T)

faith = read.table("faith_pd.tsv", header = T)
otus = read.table("observed_otus.tsv", header = T)
shannon = read.table("shannon.tsv", header = T)

alpha = inner_join(metadata, faith, by = "SampleID") %>% 
  inner_join(otus, by = "SampleID") %>% 
  inner_join(shannon, by = "SampleID")

alpha = as.data.frame(alpha)

o<-lme(fixed=observed_otus~Repro_status + Rainfall, data=alpha, 
       random= ~1|Individual)
summary(o)
Anova(o)
summary(glht(o,linfct=mcp(Repro_status="Tukey")))
r2_nakagawa(lmer(observed_otus~Repro_status + (1|Individual), data = alpha))

s<-lme(fixed=shannon~Repro_status + Rainfall, data=alpha, random= ~1|Individual)
summary(s)
Anova(s)
summary(glht(s,linfct=mcp(Repro_status="Tukey")))

oh<-lme(fixed=observed_otus~logP_avg3*Phytoprogestin + logE_avg3 + Repro_status + 
          Rainfall, data=alpha, random= ~1|Individual, na.action = na.omit)
summary(oh)
Anova(oh)
summary(glht(oh,linfct=mcp(Repro_status="Tukey")))
r2_nakagawa(lmer(observed_otus~logP_avg3 + (1|Individual), na.action = na.omit, data = alpha))
r2_nakagawa(lmer(observed_otus~logE_avg3 + (1|Individual), na.action = na.omit, data = alpha))
r2_nakagawa(lmer(observed_otus~logP_avg3*Phytoprogestin + logE_avg3 + Repro_status + 
                   Rainfall + (1|Individual), na.action = na.omit, data = alpha))

sh<-lme(fixed=shannon~logP_avg3*Phytoprogestin + logE_avg3 + Repro_status + 
          Rainfall, data=alpha, random= ~1|Individual, na.action = na.omit)
summary(sh)
Anova(sh)
summary(glht(sh,linfct=mcp(Repro_status="Tukey")))
r2_nakagawa(lmer(shannon~logP_avg3 + (1|Individual), na.action = na.omit, data = alpha))
r2_nakagawa(lmer(shannon~logE_avg3 + (1|Individual), na.action = na.omit, data = alpha))
r2_nakagawa(lmer(shannon~logP_avg3*Phytoprogestin + logE_avg3 + Repro_status + 
                   Rainfall + (1|Individual), na.action = na.omit, data = alpha))

#Alpha diversity plots----
library(ggpubr)
library(grid)

alpha$Repro_status<-factor(alpha$Repro_status,
                           levels=c("Cycling","Pregnant","Lactating"))

margr1 = grobTree(textGrob(expression(italic("marginal "*R^2*" = 0.218")), 
                           x = 0.05, y = 0.05, hjust = 0, gp = gpar(fontsize = 10)))
margr2 = grobTree(textGrob(expression(italic("marginal "*R^2*" = 0.189")), 
                           x = 0.05, y = 0.05, hjust = 0, gp = gpar(fontsize = 10)))

a = ggboxplot(alpha, x = "Repro_status", y = "shannon", fill = "Repro_status") +
  scale_fill_manual(values = c("#9fc7eb", "#d44b8e", "#8e79ba"), 
                    name = "Reproductive Status", 
                    breaks = c("Cycling", "Pregnant", 
                               "Lactating"),
                    labels = c("Cycling", "Pregnant", 
                               "Lactating")) 
a = ggpar(a, xlab = "Reproductive Status", xlegend.title = "Reproductive Status", 
          legend = "none", ylab = "Shannon Diversity Index") +
  rremove("xlab") +
  stat_compare_means(comparisons = list(c("Pregnant", "Cycling")), 
                     label = "p.signif",  symnum.args = 
                       list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.1, 1),
                            symbols = c("****", "***", "**", "*", "ns"))) 

b = ggscatter(alpha, x = "logP_avg3", y = "shannon", add = "reg.line", 
              color = "Phytoprogestin", add.params = list(color = "black")) +
  scale_color_manual(values = c("#c4e6ff", "#2d516c"))
b = ggpar(b, xlab = "Fecal Progestin (log)", ylab = "Shannon Diversity Index",
          legend.title = "Phytoprogestin\nPeriod", legend = "right") + 
  annotation_custom(margr1)


c = ggscatter(alpha, x = "logE_avg3", y = "shannon", add = "reg.line", 
              color = "Phytoprogestin", add.params = list(color = "black")) +
  scale_color_manual(values = c("#b0b0b0", "#4d4d4d"))
c = ggpar(c, xlab = "Fecal Estrogen (log)", ylab = "Shannon Diversity Index",
          legend.title = "Phytoprogestin\nPeriod", legend = "right") + 
  annotation_custom(margr1)

d = ggboxplot(alpha, x = "Repro_status", y = "observed_otus", 
              fill = "Repro_status") +
  scale_fill_manual(values = c("#9fc7eb", "#d44b8e", "#8e79ba"), 
                    name = "Reproductive Status", 
                    breaks = c("Cycling", "Pregnant", 
                               "Lactating"),
                    labels = c("Cycling", "Pregnant", 
                               "Lactating")) 
d = ggpar(d, xlab = "Reproductive Status", xlegend.title = "Reproductive Status", 
          legend = "none", ylab = "Observed OTUs") +
  rremove("xlab") +
  stat_compare_means(comparisons = list(c("Pregnant", "Cycling")), 
                     label = "p.signif",  symnum.args = 
                       list(cutpoints = c(0, 0.0001, 0.001, 0.05, 0.1, 1),
                            symbols = c("****", "***", "**", "*", "ns")))

e = ggscatter(alpha, x = "logP_avg3", y = "observed_otus", add = "reg.line", 
              color = "Phytoprogestin", add.params = list(color = "black")) +
  scale_color_manual(values = c("#c4e6ff", "#2d516c"))
e = ggpar(e, xlab = "Fecal Progestin (log)", ylab = "Observed OTUs",
          legend.title = "Phytoprogestin\nPeriod", legend = "right") + 
  annotation_custom(margr2)


f = ggscatter(alpha, x = "logE_avg3", y = "observed_otus", add = "reg.line", 
              color = "Phytoprogestin", add.params = list(color = "black")) +
  scale_color_manual(values = c("#b0b0b0", "#4d4d4d")) 
f = ggpar(f, xlab = "Fecal Estrogen (log)", ylab = "Observed OTUs",
          legend.title = "Phytoprogestin\nPeriod", legend = "right") + 
  annotation_custom(margr2)

tiff(file="alpha_rare_combined_forpub_box.tif", res=300, width=12, height=6, 
     units="in")
ggarrange(a, c, b, d, f, e, widths = c(3.5,5,5,3.5,5,5), ncol = 3, 
          nrow = 2, align = "h")
dev.off()

#PERMANOVAs----
#Import data
#sort by sampleid prior to importing

unweighted = as.dist(read.table("unweighted.tsv", header = T))
weighted = as.dist(read.table("weighted.tsv", header = T))
unweighted_bgroup = as.dist(read.table("unweighted-distance-matrix-bgroup.tsv", 
                                       header = T))
weighted_bgroup = as.dist(read.table("weighted-distance-matrix-bgroup.tsv", 
                                     header = T))

metadata_beta = read.csv("leaf_map_r_8000.csv", header = T)
metadata_beta_hormones = read.csv("leaf_map_r_8000_bgroup.csv", header = T)

#PERMANOVAs
library(vegan)

adonis(unweighted~Repro_status + Rainfall, strata=metadata_beta$Individual, 
       data=metadata_beta, permutations=5000)

adonis(weighted~Repro_status + Rainfall, strata=metadata_beta$Individual,
       data=metadata_beta, permutations=5000)

adonis(unweighted_bgroup~logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
       strata=metadata_beta_hormones$Individual, 
       data=metadata_beta_hormones, permutations=5000)

adonis(weighted_bgroup~logP_avg3 + logE_avg3 + Repro_status + Rainfall + total_per, 
       strata=metadata_beta_hormones_p$Individual, 
       data=metadata_beta_hormones, permutations=5000)

#Pairwise PERMANOVAs-----
library(pairwiseAdonis)

pairwise.adonis(unweighted, factors = metadata_beta$Repro_status, perm = 5000, 
                p.adjust.m='holm')
pairwise.adonis(weighted, factors = metadata_beta$Repro_status, perm = 5000, 
                p.adjust.m='holm')

#NMDS plots----
library(ggplot2)
library(vegan)

metadata_beta$Repro_status<-factor(metadata_beta$Repro_status,levels=
                                     c("Cycling","Pregnant","Lactating"))
metadata_beta_hormones$Repro_status<-factor(metadata_beta_hormones$Repro_status,
                                            levels=c("Cycling","Pregnant",
                                                     "Lactating"))

#Weighted full sample set
weighted.mds <- metaMDS(weighted, trymax=500)
weighted.mds.points <- weighted.mds$points
weighted.mds.points2 <- merge(x = weighted.mds.points, y = metadata_beta, 
                              by.x = "row.names", by.y = "SampleID")
grob1 = grobTree(textGrob(expression(italic(atop(R^2*" = 0.023", "p = 0.444"))), 
                                     x = 0.85, y = 0.1, hjust = 0, 
                                     gp = gpar(fontsize = 12)))
weighted_nmds <- ggplot(weighted.mds.points2, aes(x = MDS1, y = MDS2, 
                                                  color = Repro_status, 
                                                  shape = Season)) + 
  geom_point(size=4) + scale_color_manual(name = "Reproductive\nStatus", 
                                          values = c("#9fc7eb", "#d44b8e", 
                                                     "#8e79ba")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.key=element_blank()) + 
  theme(axis.title.x=element_text(size=rel(2)), 
        axis.title.y=element_text(size=rel(2)),
        plot.title = element_text(size=rel(3)),
        legend.title = element_text(size=rel(2)),
        legend.text = element_text(size = rel(1.8))) + 
  ggtitle("Weighted UniFrac") +
  stat_ellipse(aes(x = MDS1, y = MDS2, group = Repro_status, color = Repro_status), 
               type = "t") + annotation_custom(grob1)

#Unweighted full sample set
unweighted.mds <- metaMDS(unweighted, trymax=500)
unweighted.mds.points <- unweighted.mds$points
unweighted.mds.points2 <- merge(x = unweighted.mds.points, y = metadata_beta, 
                                by.x = "row.names", by.y = "SampleID")
grob2 = grobTree(textGrob(expression(italic(atop(R^2*" = 0.032", "p = 0.049"))), 
                          x = 0.85, y = 0.1, hjust = 0, 
                          gp = gpar(fontsize = 12)))
unweighted_nmds <- ggplot(unweighted.mds.points2, aes(x = MDS1, y = MDS2, 
                                                      color = Repro_status, 
                                                      shape = Season)) + 
  geom_point(size=4) + scale_color_manual(name = "Reproductive Status", 
                                          values = c("#9fc7eb", "#d44b8e", 
                                                     "#8e79ba")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.key=element_blank()) + 
  theme(axis.title.x=element_text(size=rel(2)), 
        axis.title.y=element_text(size=rel(2)),
        plot.title = element_text(size=rel(3)),
        legend.title = element_text(size=rel(2)),
        legend.text = element_text(size = rel(1.8))) + 
  ggtitle("Unweighted UniFrac") +
  stat_ellipse(aes(x = MDS1, y = MDS2, group = Repro_status, color = Repro_status), 
               type = "t") + annotation_custom(grob2)

#Weighted samples with hormones
weighted2.mds <- metaMDS(weighted_bgroup, trymax=500)
weighted2.mds.points <- weighted2.mds$points
weighted2.mds.points2 <- merge(x = weighted2.mds.points, 
                               y = metadata_beta_hormones, 
                               by.x = "row.names", by.y = "SampleID")
grob3 = grobTree(textGrob(expression(italic(atop(R^2*" = 0.043", "p = 0.010"))), 
                          x = 0.85, y = 0.1, hjust = 0, 
                          gp = gpar(fontsize = 12)))
weighted2_nmds <- ggplot(weighted2.mds.points2, aes(x = MDS1, y = logP_avg3, 
                                                    color = logP_avg3, 
                                                    shape = Season)) + 
  geom_point(size=4, color = c("#2d516c")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.key=element_blank()) + 
  theme(axis.title.x=element_text(size=rel(2)), 
        axis.title.y=element_text(size=rel(2)),
        plot.title = element_text(size=rel(3)),
        legend.title = element_text(size=rel(2)),
        legend.text = element_text(size = rel(1.8))) + 
  labs(y = "Fecal Progestin (log)") +
  ggtitle("Weighted UniFrac") + annotation_custom(grob3)

#Unweighted samples with hormones
unweighted2.mds <- metaMDS(unweighted_bgroup, trymax=500)
unweighted2.mds.points <- unweighted2.mds$points
unweighted2.mds.points2 <- merge(x = unweighted2.mds.points, 
                                 y = metadata_beta_hormones, 
                                 by.x = "row.names", by.y = "SampleID")
grob4 = grobTree(textGrob(expression(italic(atop(R^2*" = 0.029", "p = 0.004"))), 
                          x = 0.85, y = 0.1, hjust = 0, 
                          gp = gpar(fontsize = 12)))
unweighted2_nmds <- ggplot(unweighted2.mds.points2, aes(x = MDS1, y = logP_avg3, 
                                                        color = logP_avg3, 
                                                        shape = Season)) + 
  geom_point(size=4, color = c("#2d516c")) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
        legend.key=element_blank()) + 
  theme(axis.title.x=element_text(size=rel(2)), 
        axis.title.y=element_text(size=rel(2)),
        plot.title = element_text(size=rel(3)),
        legend.title = element_text(size=rel(2)),
        legend.text = element_text(size = rel(1.8))) + 
  labs(y = "Fecal Progestin (log)") +
  ggtitle("Unweighted UniFrac") + annotation_custom(grob4)

#Combined NMDS plots
library(cowplot)
legend1 = get_legend(weighted_nmds + theme(legend.box.margin = margin(0, 0, 0, 12)))
tiff(file = "combined_NMDS_plots.tif", res = 300, width=18, height=13, units="in")
row1 = plot_grid(unweighted_nmds + theme(legend.position = "none"), 
                 weighted_nmds + theme(legend.position = "none"),
                 align = "h", axis = "tb", nrow = 1, ncol = 2, labels = c('A', 'B'), 
                 label_size = 20)
row2 = plot_grid(unweighted2_nmds + theme(legend.position = "none"), 
                 weighted2_nmds + theme(legend.position = "none"),
                 align = "h", axis = "tb", nrow = 1, ncol = 2, labels = c('C', 'D'),
                 label_size = 20)
plot_grid(row1, legend1, row2, nrow = 2, ncol = 2, rel_widths = c(3, 0.5, 3))
dev.off()

#Phyla GLMMs----
library(glmmTMB)
library(multcomp)
library(car)
library(tidyverse)

glht_glmmTMB <- function (model, ..., component="cond") {
  glht(model, ...,
       coef. = function(x) fixef(x)[[component]],
       vcov. = function(x) vcov(x)[[component]],
       df = NULL)
}
modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
                               vcov. = function(x) vcov(x)[[component]],
                               df = NULL, component="cond", ...) {
  multcomp:::modelparm.default(model, coef. = coef., vcov. = vcov.,
                               df = df, ...)
}

metadata_phyla = read.csv("leaf_map_r_phyla.csv", header = T)
seqs = read.csv("phyla.csv", header = T)

phyla = inner_join(seqs, metadata_phyla, by = "SampleID") 

acti = glmmTMB(Actinobacteria ~ Repro_status + Rainfall + (1|Individual), 
               data = phyla, family = nbinom2)
summary(acti)
Anova(acti)
summary(glht(acti,linfct=mcp(Repro_status="Tukey")))

acti = glmmTMB(Actinobacteria ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + 
                 Rainfall + (1|Individual), data = phyla, family = nbinom2, 
               na.action = na.omit)
summary(acti)
Anova(acti)
summary(glht(acti,linfct=mcp(Repro_status="Tukey")))

bact = glmmTMB(Bacteroidetes ~ Repro_status + Rainfall + (1|Individual), 
               data = phyla, family = nbinom2)
summary(bact)
Anova(bact)
summary(glht(bact,linfct=mcp(Repro_status="Tukey")))

bact = glmmTMB(Bacteroidetes ~ logP_avg3*Pperiod + logE_avg3 + Repro_status +
                 Rainfall + (1|Individual), data = phyla, family = nbinom2, 
               na.action = na.omit)
summary(bact)
Anova(bact)
summary(glht(bact,linfct=mcp(Repro_status="Tukey")))

firm = glmmTMB(Firmicutes ~ Repro_status + Rainfall + (1|Individual), 
               data = phyla, family = nbinom2)
summary(firm)
Anova(firm)
summary(glht(firm,linfct=mcp(Repro_status="Tukey")))

firm = glmmTMB(Firmicutes ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + 
                 Rainfall + (1|Individual), data = phyla, family = nbinom2, 
               na.action = na.omit)
summary(firm)
Anova(firm)
summary(glht(firm,linfct=mcp(Repro_status="Tukey")))

planc = glmmTMB(Planctomycetes ~ Repro_status + Rainfall + (1|Individual), 
                data = phyla, family = nbinom2)
summary(planc)
Anova(planc)
summary(glht(planc,linfct=mcp(Repro_status="Tukey")))

planc = glmmTMB(Planctomycetes ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + 
                  Rainfall + (1|Individual), data = phyla, family = nbinom2, 
                na.action = na.omit)
summary(planc)
Anova(planc)
summary(glht(planc,linfct=mcp(Repro_status="Tukey")))

prot = glmmTMB(Proteobacteria ~ Repro_status + Rainfall + (1|Individual), 
               data = phyla, family = nbinom2)
summary(prot)
Anova(prot)
summary(glht(prot,linfct=mcp(Repro_status="Tukey")))

prot = glmmTMB(Proteobacteria ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + 
                 Rainfall + (1|Individual), data = phyla, family = nbinom2, 
               na.action = na.omit)
summary(prot)
Anova(prot)
summary(glht(prot,linfct=mcp(Repro_status="Tukey")))

spiro = glmmTMB(Spirochaetes ~ Repro_status + Rainfall + (1|Individual), 
                data = phyla, family = nbinom2)
summary(spiro)
Anova(spiro)
summary(glht(spiro,linfct=mcp(Repro_status="Tukey")))

spiro = glmmTMB(Spirochaetes ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + 
                  Rainfall + (1|Individual), data = phyla, family = nbinom2, 
                na.action = na.omit)
summary(spiro)
Anova(spiro)
summary(glht(spiro,linfct=mcp(Repro_status="Tukey")))

tene = glmmTMB(Tenericutes ~ Repro_status + Rainfall + (1|Individual), 
               data = phyla, family = nbinom2)
summary(tene)
Anova(tene)
summary(glht(tene,linfct=mcp(Repro_status="Tukey")))

tene = glmmTMB(Tenericutes ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + 
                 Rainfall + (1|Individual), data = phyla, family = nbinom2, 
               na.action = na.omit)
summary(tene)
Anova(tene)
summary(glht(tene,linfct=mcp(Repro_status="Tukey")))

#Family GLMMs----
library(glmmTMB)
library(multcomp)
library(car)
library(tidyverse)

glht_glmmTMB <- function (model, ..., component="cond") {
  glht(model, ...,
       coef. = function(x) fixef(x)[[component]],
       vcov. = function(x) vcov(x)[[component]],
       df = NULL)
}
modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
                               vcov. = function(x) vcov(x)[[component]],
                               df = NULL, component="cond", ...) {
  multcomp:::modelparm.default(model, coef. = coef., vcov. = vcov.,
                               df = df, ...)
}

family_relab = read.table("feature-table-level5.tsv", header = T)
metadata = read.csv("leaf_map_r_phyla.csv", header = T)

#Without hormones

family = dplyr::select(metadata, "SampleID":"Pperiod") %>% 
  inner_join(family_relab, by = "SampleID")

repro_estimatematrix = mat.or.vec(128,2)
repro_pvaluematrix = mat.or.vec(128,2)
rainfall_estimatematrix = mat.or.vec(128,2)
rainfall_pvaluematrix = mat.or.vec(128,2)

for(i in 11:137) {
  variable = family[,i]
  b<-try(glmmTMB(variable ~ Repro_status + Rainfall + (1|Individual), 
                 data = family, family = nbinom2))
  anova = Anova(b)
  repro_estimatematrix[i-10,2] = anova[1,1]
  repro_pvaluematrix[i-10,2] = anova[1,3]
  rainfall_estimatematrix[i-10,2] = anova[2,1]
  rainfall_pvaluematrix[i-10,2] = anova[2,3]
  repro_estimatematrix[i-10,1] = names(family)[i]
  repro_pvaluematrix[i-10,1] = names(family)[i]
  rainfall_estimatematrix[i-10,1] = names(family)[i]
  rainfall_pvaluematrix[i-10,1] = names(family)[i]
}

family_repro = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                         as.data.frame(repro_pvaluematrix[,2]))
write.csv(family_repro, "family_repro.csv")
family_repro_nona = filter(family_repro, family_repro[,3] != "NaN")
family_repro_nona[,3] = as.numeric(as.character(family_repro_nona[,3]))
family_repro_corrected = bind_cols(family_repro_nona, 
                                   as.data.frame(fdrtool(family_repro_nona[,3], 
                                                         statistic = "pvalue", 
                                                         plot = F)))
write.csv(family_repro_corrected, "family_repro_corrected.csv")

family_rainfall = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                            as.data.frame(rainfall_pvaluematrix[,2]))
write.csv(family_rainfall, "family_rainfall.csv")
family_rainfall_nona = filter(family_rainfall, family_rainfall[,3] != "NaN")
family_rainfall_nona[,3] = as.numeric(as.character(family_rainfall_nona[,3]))
family_rainfall_corrected = bind_cols(family_rainfall_nona, 
                                      as.data.frame(fdrtool(family_rainfall_nona[,3], 
                                                            statistic = "pvalue", 
                                                            plot = F)))
write.csv(family_rainfall_corrected, "family_rainfall_corrected.csv")

#With hormones

family_hormones = dplyr::select(metadata, "SampleID":"Pperiod", "logP_avg3", 
                                "logE_avg3") %>% filter(logP_avg3 > 0) %>% 
  inner_join(family_relab, by = "SampleID") 

repro_estimatematrix = mat.or.vec(128,2)
repro_pvaluematrix = mat.or.vec(128,2)
rainfall_estimatematrix = mat.or.vec(128,2)
rainfall_pvaluematrix = mat.or.vec(128,2)
p_estimatematrix = mat.or.vec(128,2)
p_pvaluematrix = mat.or.vec(128,2)
e_estimatematrix = mat.or.vec(128,2)
e_pvaluematrix = mat.or.vec(128,2)
ebypp_estimatematrix = mat.or.vec(128,2)
ebypp_pvaluematrix = mat.or.vec(128,2)
pp_estimatematrix = mat.or.vec(128,2)
pp_pvaluematrix = mat.or.vec(128,2)

for(i in 140:140) {
  variable = family_hormones[,i]
  b = try(glmmTMB(variable ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + 
                    Rainfall + (1|Individual), data = family_hormones, 
                  family = nbinom2))
  anova = Anova(b)
  repro_estimatematrix[i-13,2] = anova[4,1]
  repro_pvaluematrix[i-13,2] = anova[4,3]
  rainfall_estimatematrix[i-13,2] = anova[5,1]
  rainfall_pvaluematrix[i-13,2] = anova[5,3]
  p_estimatematrix[i-13,2] = anova[1,1]
  p_pvaluematrix[i-13,2] = anova[1,3]
  e_estimatematrix[i-13,2] = anova[3,1]
  e_pvaluematrix[i-13,2] = anova[3,3]
  pp_estimatematrix[i-13,2] = anova[2,1]
  pp_pvaluematrix[i-13,2] = anova[2,3]
  ebypp_estimatematrix[i-13,2] = anova[6,1]
  ebypp_pvaluematrix[i-13,2] = anova[6,3]
  repro_estimatematrix[i-13,1] = names(family)[i]
  repro_pvaluematrix[i-13,1] = names(family)[i]
  rainfall_estimatematrix[i-13,1] = names(family)[i]
  rainfall_pvaluematrix[i-13,1] = names(family)[i]
  p_estimatematrix[i-13,1] = names(family)[i]
  p_pvaluematrix[i-13,1] = names(family)[i]
  e_estimatematrix[i-13,1] = names(family)[i]
  e_pvaluematrix[i-13,1] = names(family)[i]
  pp_estimatematrix[i-13,1] = names(family)[i]
  pp_pvaluematrix[i-13,1] = names(family)[i]
  ebypp_estimatematrix[i-13,1] = names(family)[i]
  ebypp_pvaluematrix[i-13,1] = names(family)[i]
}

family_repro_hormones = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                                  as.data.frame(repro_pvaluematrix[,2]))
write.csv(family_repro_hormones, "family_repro_hormones.csv")
family_repro_hormones_nona = filter(family_repro_hormones, 
                                    family_repro_hormones[,3] != "NaN")
family_repro_hormones_nona[,3] = as.numeric(as.character(family_repro_hormones_nona[,3]))
family_repro_hormones_corrected = bind_cols(family_repro_hormones_nona, 
                                            as.data.frame(fdrtool(family_repro_hormones_nona[,3], 
                                                                  statistic = "pvalue", plot = F)))
write.csv(family_repro_hormones_corrected, "family_repro_hormones_corrected.csv")

family_rainfall_hormones = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                                     as.data.frame(rainfall_pvaluematrix[,2]))
write.csv(family_rainfall_hormones, "family_rainfall_hormones.csv")
family_rainfall_hormones_nona = filter(family_rainfall_hormones, 
                                       family_rainfall_hormones[,3] != "NaN")
family_rainfall_hormones_nona[,3] = as.numeric(as.character(family_rainfall_hormones_nona[,3]))
family_rainfall_hormones_corrected = bind_cols(family_rainfall_hormones_nona, 
                                               as.data.frame(fdrtool(family_rainfall_hormones_nona[,3], 
                                                                     statistic = "pvalue", plot = F)))
write.csv(family_rainfall_hormones_corrected, "family_rainfall_hormones_corrected.csv")

family_prog_hormones = bind_cols(as.data.frame(p_estimatematrix[,1:2]), 
                                 as.data.frame(p_pvaluematrix[,2]))
write.csv(family_prog_hormones, "family_prog_hormones.csv")
family_prog_hormones_nona = filter(family_prog_hormones, 
                                   family_prog_hormones[,3] != "NaN")
family_prog_hormones_nona[,3] = as.numeric(as.character(family_prog_hormones_nona[,3]))
family_prog_hormones_corrected = bind_cols(family_prog_hormones_nona, 
                                           as.data.frame(fdrtool(family_prog_hormones_nona[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(family_prog_hormones_corrected, "family_prog_hormones_corrected.csv")

family_estr_hormones = bind_cols(as.data.frame(e_estimatematrix[,1:2]), 
                                 as.data.frame(e_pvaluematrix[,2]))
write.csv(family_estr_hormones, "family_estr_hormones.csv")
family_estr_hormones_nona = filter(family_estr_hormones, 
                                   family_estr_hormones[,3] != "NaN")
family_estr_hormones_nona[,3] = as.numeric(as.character(family_estr_hormones_nona[,3]))
family_estr_hormones_corrected = bind_cols(family_estr_hormones_nona, 
                                           as.data.frame(fdrtool(family_estr_hormones_nona[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(family_estr_hormones_corrected, "family_estr_hormones_corrected.csv")

family_pperiod_hormones = bind_cols(as.data.frame(pp_estimatematrix[,1:2]), 
                                    as.data.frame(pp_pvaluematrix[,2]))
write.csv(family_pperiod_hormones, "family_pperiod_hormones.csv")
family_pperiod_hormones_nona = filter(family_pperiod_hormones, 
                                      family_period_hormones[,3] != "NaN")
family_pperiod_hormones_nona[,3] = as.numeric(as.character(family_pperiod_hormones_nona[,3]))
family_pperiod_hormones_corrected = bind_cols(family_period_hormones_nona, 
                                              as.data.frame(fdrtool(family_pperiod_hormones_nona[,3], 
                                                                    statistic = "pvalue", plot = F)))
write.csv(family_pperiod_hormones_corrected, "family_pperiod_hormones_corrected.csv")

family_ebypp_hormones = bind_cols(as.data.frame(ebypp_estimatematrix[,1:2]), 
                                  as.data.frame(ebypp_pvaluematrix[,2]))
write.csv(family_ebypp_hormones, "family_ebypp_hormones.csv")
family_ebypp_hormones_nona = filter(family_ebypp_hormones, 
                                    family_ebypp_hormones[,3] != "NaN")
family_ebypp_hormones_nona[,3] = as.numeric(as.character(family_ebypp_hormones_nona[,3]))
family_ebypp_hormones_corrected = bind_cols(family_ebypp_hormones_nona, 
                                            as.data.frame(fdrtool(family_ebypp_hormones_nona[,3], 
                                                                  statistic = "pvalue", plot = F)))
write.csv(family_ebypp_hormones_corrected, "family_ebypp_hormones_corrected.csv")

#Genus GLMMs----
library(glmmTMB)
library(multcomp)
library(car)
library(tidyverse)

glht_glmmTMB <- function (model, ..., component="cond") {
  glht(model, ...,
       coef. = function(x) fixef(x)[[component]],
       vcov. = function(x) vcov(x)[[component]],
       df = NULL)
}
modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
                               vcov. = function(x) vcov(x)[[component]],
                               df = NULL, component="cond", ...) {
  multcomp:::modelparm.default(model, coef. = coef., vcov. = vcov.,
                               df = df, ...)
}

genus_relab = read.csv("feature-table-level6.csv", header = T)
metadata = read.csv("leaf_map_r_phyla.csv", header = T)

#Without hormones

genus = dplyr::select(metadata, 1:10) %>%inner_join(genus_relab, by = "SampleID")

repro_estimatematrix = mat.or.vec(185,2)
repro_pvaluematrix = mat.or.vec(185,2)
rainfall_estimatematrix = mat.or.vec(185,2)
rainfall_pvaluematrix = mat.or.vec(185,2)

for(i in 106:194) {
  variable = genus[,i]
  b<-try(glmmTMB(variable ~ Repro_status + Rainfall + (1|Individual), 
                 data = genus, family = nbinom2))
  anova = Anova(b)
  repro_estimatematrix[i-10,2] = anova[1,1]
  repro_pvaluematrix[i-10,2] = anova[1,3]
  rainfall_estimatematrix[i-10,2] = anova[2,1]
  rainfall_pvaluematrix[i-10,2] = anova[2,3]
  repro_estimatematrix[i-10,1] = names(genus)[i]
  repro_pvaluematrix[i-10,1] = names(genus)[i]
  rainfall_estimatematrix[i-10,1] = names(genus)[i]
  rainfall_pvaluematrix[i-10,1] = names(genus)[i]
}

genus_repro = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                        as.data.frame(repro_pvaluematrix[,2]))
write.csv(genus_repro, "genus_repro.csv")
genus_repro_nona = filter(genus_repro, genus_repro[,3] != "NaN")
genus_repro_nona[,3] = as.numeric(as.character(genus_repro_nona[,3]))
genus_repro_corrected = bind_cols(genus_repro_nona, 
                                  as.data.frame(fdrtool(genus_repro_nona[,3], 
                                                        statistic = "pvalue", plot = F)))
write.csv(genus_repro_corrected, "genus_repro_corrected.csv")

genus_rainfall = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                           as.data.frame(rainfall_pvaluematrix[,2]))
write.csv(genus_rainfall, "genus_rainfall.csv")
genus_rainfall_nona = filter(genus_rainfall, genus_rainfall[,3] != "NaN")
genus_rainfall_nona[,3] = as.numeric(as.character(genus_rainfall_nona[,3]))
genus_rainfall_corrected = bind_cols(genus_rainfall_nona, 
                                     as.data.frame(fdrtool(genus_rainfall_nona[,3], 
                                                           statistic = "pvalue", plot = F)))
write.csv(genus_rainfall_corrected, "genus_rainfall_corrected.csv")

#With hormones

genus_hormones = dplyr::select(metadata, "SampleID":"Rainfall", "logP_avg3", 
                               "logE_avg3") %>% filter(logP_avg3 > 0) %>% 
  inner_join(genus_relab, by = "SampleID") 

repro_estimatematrix = mat.or.vec(185,2)
repro_pvaluematrix = mat.or.vec(185,2)
rainfall_estimatematrix = mat.or.vec(185,2)
rainfall_pvaluematrix = mat.or.vec(185,2)
p_estimatematrix = mat.or.vec(185,2)
p_pvaluematrix = mat.or.vec(185,2)
e_estimatematrix = mat.or.vec(185,2)
e_pvaluematrix = mat.or.vec(185,2)

for(i in 13:196) {
  variable = genus_hormones[,i]
  b = try(glmmTMB(variable ~ logP_avg3 + logE_avg3 + Repro_status + Rainfall + 
                    (1|Individual), data = genus_hormones, family = nbinom2))
  anova = Anova(b)
  repro_estimatematrix[i-12,2] = anova[3,1]
  repro_pvaluematrix[i-12,2] = anova[3,3]
  rainfall_estimatematrix[i-12,2] = anova[4,1]
  rainfall_pvaluematrix[i-12,2] = anova[4,3]
  p_estimatematrix[i-12,2] = anova[1,1]
  p_pvaluematrix[i-12,2] = anova[1,3]
  e_estimatematrix[i-12,2] = anova[2,1]
  e_pvaluematrix[i-12,2] = anova[2,3]
  repro_estimatematrix[i-12,1] = names(genus)[i]
  repro_pvaluematrix[i-12,1] = names(genus)[i]
  rainfall_estimatematrix[i-12,1] = names(genus)[i]
  rainfall_pvaluematrix[i-12,1] = names(genus)[i]
  p_estimatematrix[i-12,1] = names(genus)[i]
  p_pvaluematrix[i-12,1] = names(genus)[i]
  e_estimatematrix[i-12,1] = names(genus)[i]
  e_pvaluematrix[i-12,1] = names(genus)[i]
}

genus_repro_hormones = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                                 as.data.frame(repro_pvaluematrix[,2]))
write.csv(genus_repro_hormones, "genus_repro_hormones.csv")
genus_repro_hormones_nona = filter(genus_repro_hormones, 
                                   genus_repro_hormones[,3] != "NaN")
genus_repro_hormones_nona[,3] = as.numeric(as.character(genus_repro_hormones_nona[,3]))
genus_repro_hormones_corrected = bind_cols(genus_repro_hormones_nona, 
                                           as.data.frame(fdrtool(genus_repro_hormones_nona[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(genus_repro_hormones_corrected, "genus_repro_hormones_corrected.csv")

genus_rainfall_hormones = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                                    as.data.frame(rainfall_pvaluematrix[,2]))
write.csv(genus_rainfall_hormones, "genus_rainfall_hormones.csv")
genus_rainfall_hormones_nona = filter(genus_rainfall_hormones, 
                                      genus_rainfall_hormones[,3] != "NaN")
genus_rainfall_hormones_nona[,3] = as.numeric(as.character(genus_rainfall_hormones_nona[,3]))
genus_rainfall_hormones_corrected = bind_cols(genus_rainfall_hormones_nona, 
                                              as.data.frame(fdrtool(genus_rainfall_hormones_nona[,3], 
                                                                    statistic = "pvalue", plot = F)))
write.csv(genus_rainfall_hormones_corrected, "genus_rainfall_hormones_corrected.csv")

genus_prog_hormones = bind_cols(as.data.frame(p_estimatematrix[,1:2]), 
                                as.data.frame(p_pvaluematrix[,2]))
write.csv(genus_prog_hormones, "genus_prog_hormones.csv")
genus_prog_hormones_nona = filter(genus_prog_hormones, 
                                  genus_prog_hormones[,3] != "NaN")
genus_prog_hormones_nona[,3] = as.numeric(as.character(genus_prog_hormones_nona[,3]))
genus_prog_hormones_corrected = bind_cols(genus_prog_hormones_nona, 
                                          as.data.frame(fdrtool(genus_prog_hormones_nona[,3], 
                                                                statistic = "pvalue", plot = F)))
write.csv(genus_prog_hormones_corrected, "genus_prog_hormones_corrected.csv")

genus_estr_hormones = bind_cols(as.data.frame(e_estimatematrix[,1:2]), 
                                as.data.frame(e_pvaluematrix[,2]))
write.csv(genus_estr_hormones, "genus_estr_hormones.csv")
genus_estr_hormones_nona = filter(genus_estr_hormones, 
                                  genus_estr_hormones[,3] != "NaN")
genus_estr_hormones_nona[,3] = as.numeric(as.character(genus_estr_hormones_nona[,3]))
genus_estr_hormones_corrected = bind_cols(genus_estr_hormones_nona, 
                                          as.data.frame(fdrtool(genus_estr_hormones_nona[,3], 
                                                                statistic = "pvalue", plot = F)))
write.csv(genus_estr_hormones_corrected, "genus_estr_hormones_corrected.csv")


#Phyla presence-absence GLMMs----
library(multcomp)
library(lme4)
library(car)
library(tidyverse)

metadata_phyla = read.csv("leaf_map_r_phyla.csv", header = T)
seqs = read.csv("phyla_pa.csv", header = T)

phyla = inner_join(seqs, metadata_phyla, by = "SampleID") 

acti = glmer(Actinobacteria ~ Repro_status + Rainfall + (1|Individual), 
             data = phyla, family = binomial)
summary(acti)
Anova(acti)
summary(glht(acti,linfct=mcp(Repro_status="Tukey")))

acti = glmer(Actinobacteria ~ logP_avg3 + logE_avg3 + Repro_status + Rainfall + (1|Individual), 
             data = phyla, family = binomial, na.action = na.omit)
summary(acti)
Anova(acti)
summary(glht(acti,linfct=mcp(Repro_status="Tukey")))

bact = glmer(Bacteroidetes ~ Repro_status + Rainfall + (1|Individual), 
             data = phyla, family = binomial)
summary(bact)
Anova(bact)
summary(glht(bact,linfct=mcp(Repro_status="Tukey")))

bact = glmer(Bacteroidetes ~ logP_avg3 + logE_avg3 + Repro_status + Rainfall + (1|Individual), 
             data = phyla, family = binomial, na.action = na.omit)
summary(bact)
Anova(bact)
summary(glht(bact,linfct=mcp(Repro_status="Tukey")))

firm = glmer(Firmicutes ~ Repro_status + Rainfall + (1|Individual), 
             data = phyla, family = binomial)
summary(firm)
Anova(firm)
summary(glht(firm,linfct=mcp(Repro_status="Tukey")))

firm = glmer(Firmicutes ~ logP_avg3 + logE_avg3 + Repro_status + Rainfall + (1|Individual), 
             data = phyla, family = binomial, na.action = na.omit)
summary(firm)
Anova(firm)
summary(glht(firm,linfct=mcp(Repro_status="Tukey")))

planc = glmer(Planctomycetes ~ Repro_status + Rainfall + (1|Individual), 
              data = phyla, family = binomial)
summary(planc)
Anova(planc)
summary(glht(planc,linfct=mcp(Repro_status="Tukey")))

planc = glmer(Planctomycetes ~ logP_avg3 + logE_avg3 + Repro_status + Rainfall + (1|Individual), 
              data = phyla, family = binomial, na.action = na.omit)
summary(planc)
Anova(planc)
summary(glht(planc,linfct=mcp(Repro_status="Tukey")))

prot = glmer(Proteobacteria ~ Repro_status + Rainfall + (1|Individual), 
             data = phyla, family = binomial)
summary(prot)
Anova(prot)
summary(glht(prot,linfct=mcp(Repro_status="Tukey")))

prot = glmer(Proteobacteria ~ logP_avg3 + logE_avg3 + Repro_status + Rainfall + (1|Individual), 
             data = phyla, family = binomial, na.action = na.omit)
summary(prot)
Anova(prot)
summary(glht(prot,linfct=mcp(Repro_status="Tukey")))

spiro = glmer(Spirochaetes ~ Repro_status + Rainfall + (1|Individual), 
              data = phyla, family = binomial(link = "cloglog"))
summary(spiro)
Anova(spiro)
summary(glht(spiro,linfct=mcp(Repro_status="Tukey")))

spiro = glmer(Spirochaetes ~ logP_avg3 + logE_avg3 + Repro_status + Rainfall + (1|Individual), 
              data = phyla, family = binomial(link = "cloglog"), na.action = na.omit)
summary(spiro)
Anova(spiro)
summary(glht(spiro,linfct=mcp(Repro_status="Tukey")))

tene = glmer(Tenericutes ~ Repro_status + Rainfall + (1|Individual), 
             data = phyla, family = binomial)
summary(tene)
Anova(tene)
summary(glht(tene,linfct=mcp(Repro_status="Tukey")))

tene = glmer(Tenericutes ~ logP_avg3 + logE_avg3 + Repro_status + Rainfall + (1|Individual), 
             data = phyla, family = binomial, na.action = na.omit)
summary(tene)
Anova(tene)
summary(glht(tene,linfct=mcp(Repro_status="Tukey")))

#Family presence-absence GLMMs----
library(multcomp)
library(lme4)
library(car)
library(tidyverse)
library(fdrtool)

family_pa = read.csv("family_pa.csv", header = T)
metadata = read.csv("leaf_map_r_phyla.csv", header = T)

#Without hormones

family = dplyr::select(metadata, "SampleID":"Rainfall") %>% inner_join(family_pa, by = "SampleID")

repro_estimatematrix = mat.or.vec(128,2)
repro_pvaluematrix = mat.or.vec(128,2)
rainfall_estimatematrix = mat.or.vec(128,2)
rainfall_pvaluematrix = mat.or.vec(128,2)

for(i in 136:137) {
  variable = family[,i]
  b<-glmer(variable ~ Repro_status + Rainfall + (1|Individual), 
           data = family, family = binomial)
  anova = Anova(b)
  repro_estimatematrix[i-10,2] = anova[1,1]
  repro_pvaluematrix[i-10,2] = anova[1,3]
  rainfall_estimatematrix[i-10,2] = anova[2,1]
  rainfall_pvaluematrix[i-10,2] = anova[2,3]
  repro_estimatematrix[i-10,1] = names(family)[i]
  repro_pvaluematrix[i-10,1] = names(family)[i]
  rainfall_estimatematrix[i-10,1] = names(family)[i]
  rainfall_pvaluematrix[i-10,1] = names(family)[i]
}

family_repro = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                         as.data.frame(repro_pvaluematrix[,2]))
write.csv(family_repro, "family_repro_pa.csv")
family_repro_nona = filter(family_repro, family_repro[,3] != 0)
family_repro_nona[,3] = as.numeric(as.character(family_repro_nona[,3]))
family_repro_corrected = bind_cols(family_repro_nona, 
                                   as.data.frame(fdrtool(family_repro_nona[,3], 
                                                         statistic = "pvalue", plot = F)))
write.csv(family_repro_corrected, "family_repro_pa_corrected.csv")

family_rainfall = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                            as.data.frame(rainfall_pvaluematrix[,2]))
write.csv(family_rainfall, "family_rainfall_pa.csv")
family_rainfall_nona = filter(family_rainfall, family_rainfall[,3] != 0)
family_rainfall_nona[,3] = as.numeric(as.character(family_rainfall_nona[,3]))
family_rainfall_corrected = bind_cols(family_rainfall_nona, 
                                      as.data.frame(fdrtool(family_rainfall_nona[,3], 
                                                            statistic = "pvalue", plot = F)))
write.csv(family_rainfall_corrected, "family_rainfall_pa_corrected.csv")

#With hormones

family_hormones = dplyr::select(metadata, "SampleID":"Pperiod", "logP_avg3", 
                                "logE_avg3") %>% filter(logP_avg3 > 0) %>% 
  inner_join(family_pa, by = "SampleID") 

repro_estimatematrix = mat.or.vec(128,2)
repro_pvaluematrix = mat.or.vec(128,2)
rainfall_estimatematrix = mat.or.vec(128,2)
rainfall_pvaluematrix = mat.or.vec(128,2)
p_estimatematrix = mat.or.vec(128,2)
p_pvaluematrix = mat.or.vec(128,2)
e_estimatematrix = mat.or.vec(128,2)
e_pvaluematrix = mat.or.vec(128,2)
ebypp_estimatematrix = mat.or.vec(128,2)
ebypp_pvaluematrix = mat.or.vec(128,2)
pp_estimatematrix = mat.or.vec(128,2)
pp_pvaluematrix = mat.or.vec(128,2)

for(i in 140:140) {
  variable = family_hormones[,i]
  b = glmer(variable ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + 
              Rainfall + (1|Individual), data = family_hormones, 
            family = binomial)
  anova = Anova(b)
  repro_estimatematrix[i-13,2] = anova[4,1]
  repro_pvaluematrix[i-13,2] = anova[4,3]
  rainfall_estimatematrix[i-13,2] = anova[5,1]
  rainfall_pvaluematrix[i-13,2] = anova[5,3]
  p_estimatematrix[i-13,2] = anova[1,1]
  p_pvaluematrix[i-13,2] = anova[1,3]
  e_estimatematrix[i-13,2] = anova[3,1]
  e_pvaluematrix[i-13,2] = anova[3,3]
  pp_estimatematrix[i-13,2] = anova[2,1]
  pp_pvaluematrix[i-13,2] = anova[2,3]
  ebypp_estimatematrix[i-13,2] = anova[6,1]
  ebypp_pvaluematrix[i-13,2] = anova[6,3]
  repro_estimatematrix[i-13,1] = names(family)[i]
  repro_pvaluematrix[i-13,1] = names(family)[i]
  rainfall_estimatematrix[i-13,1] = names(family)[i]
  rainfall_pvaluematrix[i-13,1] = names(family)[i]
  p_estimatematrix[i-13,1] = names(family)[i]
  p_pvaluematrix[i-13,1] = names(family)[i]
  e_estimatematrix[i-13,1] = names(family)[i]
  e_pvaluematrix[i-13,1] = names(family)[i]
  pp_estimatematrix[i-13,1] = names(family)[i]
  pp_pvaluematrix[i-13,1] = names(family)[i]
  ebypp_estimatematrix[i-13,1] = names(family)[i]
  ebypp_pvaluematrix[i-13,1] = names(family)[i]
}

family_repro_hormones = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                                  as.data.frame(repro_pvaluematrix[,2]))
write.csv(family_repro_hormones, "family_repro_hormones_pa.csv")
family_repro_hormones_nona = filter(family_repro_hormones, 
                                    family_repro_hormones[,3] != "NaN")
family_repro_hormones_nona[,3] = as.numeric(as.character(family_repro_hormones_nona[,3]))
family_repro_hormones_corrected = bind_cols(family_repro_hormones_nona, 
                                            as.data.frame(fdrtool(family_repro_hormones_nona[,3], 
                                                                  statistic = "pvalue", plot = F)))
write.csv(family_repro_hormones_corrected, "family_repro_hormones_pa_corrected.csv")

family_rainfall_hormones = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                                     as.data.frame(rainfall_pvaluematrix[,2]))
write.csv(family_rainfall_hormones, "family_rainfall_hormones_pa.csv")
family_rainfall_hormones_nona = filter(family_rainfall_hormones, 
                                       family_rainfall_hormones[,3] != "NaN")
family_rainfall_hormones_nona[,3] = as.numeric(as.character(family_rainfall_hormones_nona[,3]))
family_rainfall_hormones_corrected = bind_cols(family_rainfall_hormones_nona, 
                                               as.data.frame(fdrtool(family_rainfall_hormones_nona[,3], 
                                                                     statistic = "pvalue", plot = F)))
write.csv(family_rainfall_hormones_corrected, "family_rainfall_hormones_pa_corrected.csv")

family_prog_hormones = bind_cols(as.data.frame(p_estimatematrix[,1:2]), 
                                 as.data.frame(p_pvaluematrix[,2]))
write.csv(family_prog_hormones, "family_prog_hormones_pa.csv")
family_prog_hormones_nona = filter(family_prog_hormones, 
                                   family_prog_hormones[,3] != "NaN")
family_prog_hormones_nona[,3] = as.numeric(as.character(family_prog_hormones_nona[,3]))
family_prog_hormones_corrected = bind_cols(family_prog_hormones_nona, 
                                           as.data.frame(fdrtool(family_prog_hormones_nona[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(family_prog_hormones_corrected, "family_prog_hormones_pa_corrected.csv")

family_estr_hormones = bind_cols(as.data.frame(e_estimatematrix[,1:2]), 
                                 as.data.frame(e_pvaluematrix[,2]))
write.csv(family_estr_hormones, "family_estr_hormones_pa.csv")
family_estr_hormones_nona = filter(family_estr_hormones, 
                                   family_estr_hormones[,3] != "NaN")
family_estr_hormones_nona[,3] = as.numeric(as.character(family_estr_hormones_nona[,3]))
family_estr_hormones_corrected = bind_cols(family_estr_hormones_nona, 
                                           as.data.frame(fdrtool(family_estr_hormones_nona[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(family_estr_hormones_corrected, "family_estr_hormones_pa_corrected.csv")

family_pperiod_hormones = bind_cols(as.data.frame(pp_estimatematrix[,1:2]), 
                                    as.data.frame(pp_pvaluematrix[,2]))
write.csv(family_pperiod_hormones, "family_pperiod_hormones_pa.csv")
family_pperiod_hormones_nona = filter(family_pperiod_hormones, 
                                      family_pperiod_hormones[,3] != "NaN")
family_pperiod_hormones_nona[,3] = as.numeric(as.character(family_pperiod_hormones_nona[,3]))
family_pperiod_hormones_corrected = bind_cols(family_pperiod_hormones_nona, 
                                              as.data.frame(fdrtool(family_pperiod_hormones_nona[,3], 
                                                                    statistic = "pvalue", plot = F)))
write.csv(family_pperiod_hormones_corrected, "family_pperiod_hormones_pa_corrected.csv")

family_ebypp_hormones = bind_cols(as.data.frame(ebypp_estimatematrix[,1:2]), 
                                  as.data.frame(ebypp_pvaluematrix[,2]))
write.csv(family_ebypp_hormones, "family_ebypp_hormones_pa.csv")
family_ebypp_hormones_nona = filter(family_ebypp_hormones, 
                                    family_ebypp_hormones[,3] != "NaN")
family_ebypp_hormones_nona[,3] = as.numeric(as.character(family_ebypp_hormones_nona[,3]))
family_ebypp_hormones_corrected = bind_cols(family_ebypp_hormones_nona, 
                                            as.data.frame(fdrtool(family_ebypp_hormones_nona[,3], 
                                                                  statistic = "pvalue", plot = F)))
write.csv(family_ebypp_hormones_corrected, "family_ebypp_hormones_pa_corrected.csv")


#Genus presence-absence GLMMs----
library(multcomp)
library(lme4)
library(car)
library(tidyverse)
library(fdrtool)
genus_pa = read.csv("genus_pa.csv", header = T)
metadata = read.csv("leaf_map_r_phyla.csv", header = T)

#Without hormones

genus = dplyr::select(metadata, "SampleID":"Rainfall") %>% inner_join(genus_pa, by = "SampleID")

repro_estimatematrix = mat.or.vec(185,2)
repro_pvaluematrix = mat.or.vec(185,2)
rainfall_estimatematrix = mat.or.vec(185,2)
rainfall_pvaluematrix = mat.or.vec(185,2)

for(i in 194:194) {
  variable = genus[,i]
  b<-glmer(variable ~ Repro_status + Rainfall + (1|Individual), 
           data = genus, family = binomial)
  anova = Anova(b)
  repro_estimatematrix[i-10,2] = anova[1,1]
  repro_pvaluematrix[i-10,2] = anova[1,3]
  rainfall_estimatematrix[i-10,2] = anova[2,1]
  rainfall_pvaluematrix[i-10,2] = anova[2,3]
  repro_estimatematrix[i-10,1] = names(genus)[i]
  repro_pvaluematrix[i-10,1] = names(genus)[i]
  rainfall_estimatematrix[i-10,1] = names(genus)[i]
  rainfall_pvaluematrix[i-10,1] = names(genus)[i]
}

genus_repro = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                        as.data.frame(repro_pvaluematrix[,2]))
write.csv(genus_repro, "genus_repro_pa.csv")
genus_repro_nona = filter(genus_repro, genus_repro[,3] != 0)
genus_repro_nona[,3] = as.numeric(as.character(genus_repro_nona[,3]))
genus_repro_corrected = bind_cols(genus_repro_nona, 
                                  as.data.frame(fdrtool(genus_repro_nona[,3], 
                                                        statistic = "pvalue", plot = F)))
write.csv(genus_repro_corrected, "genus_repro_pa_corrected.csv")

genus_rainfall = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                           as.data.frame(rainfall_pvaluematrix[,2]))
write.csv(genus_rainfall, "genus_rainfall_pa.csv")
genus_rainfall_nona = filter(genus_rainfall, genus_rainfall[,3] != 0)
genus_rainfall_nona[,3] = as.numeric(as.character(genus_rainfall_nona[,3]))
genus_rainfall_corrected = bind_cols(genus_rainfall_nona, 
                                     as.data.frame(fdrtool(genus_rainfall_nona[,3], 
                                                           statistic = "pvalue", plot = F)))
write.csv(genus_rainfall_corrected, "genus_rainfall_pa_corrected.csv")

#With hormones

genus_hormones = dplyr::select(metadata, "SampleID":"Pperiod", "logP_avg3", 
                               "logE_avg3") %>% filter(logP_avg3 > 0) %>% 
  inner_join(genus_pa, by = "SampleID") 

repro_estimatematrix = mat.or.vec(185,2)
repro_pvaluematrix = mat.or.vec(185,2)
rainfall_estimatematrix = mat.or.vec(185,2)
rainfall_pvaluematrix = mat.or.vec(185,2)
p_estimatematrix = mat.or.vec(185,2)
p_pvaluematrix = mat.or.vec(185,2)
e_estimatematrix = mat.or.vec(185,2)
e_pvaluematrix = mat.or.vec(185,2)
ebypp_estimatematrix = mat.or.vec(185,2)
ebypp_pvaluematrix = mat.or.vec(185,2)
pp_estimatematrix = mat.or.vec(185,2)
pp_pvaluematrix = mat.or.vec(185,2)

for(i in 197:197) {
  variable = genus_hormones[,i]
  b = glmer(variable ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall + 
              (1|Individual), data = genus_hormones, family = binomial)
  anova = Anova(b)
  repro_estimatematrix[i-13,2] = anova[4,1]
  repro_pvaluematrix[i-13,2] = anova[4,3]
  rainfall_estimatematrix[i-13,2] = anova[5,1]
  rainfall_pvaluematrix[i-13,2] = anova[5,3]
  p_estimatematrix[i-13,2] = anova[1,1]
  p_pvaluematrix[i-13,2] = anova[1,3]
  e_estimatematrix[i-13,2] = anova[3,1]
  e_pvaluematrix[i-13,2] = anova[3,3]
  pp_estimatematrix[i-13,2] = anova[2,1]
  pp_pvaluematrix[i-13,2] = anova[2,3]
  ebypp_estimatematrix[i-13,2] = anova[6,1]
  ebypp_pvaluematrix[i-13,2] = anova[6,3]
  repro_estimatematrix[i-13,1] = names(genus_hormones)[i]
  repro_pvaluematrix[i-13,1] = names(genus_hormones)[i]
  rainfall_estimatematrix[i-13,1] = names(genus_hormones)[i]
  rainfall_pvaluematrix[i-13,1] = names(genus_hormones)[i]
  p_estimatematrix[i-13,1] = names(genus_hormones)[i]
  p_pvaluematrix[i-13,1] = names(genus_hormones)[i]
  e_estimatematrix[i-13,1] = names(genus_hormones)[i]
  e_pvaluematrix[i-13,1] = names(genus_hormones)[i]
}

genus_repro_hormones = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                                 as.data.frame(repro_pvaluematrix[,2]))
write.csv(genus_repro_hormones, "genus_repro_hormones_pa.csv")
genus_repro_hormones_nona = filter(genus_repro_hormones, 
                                   genus_repro_hormones[,3] != 0)
genus_repro_hormones_nona[,3] = as.numeric(as.character(genus_repro_hormones_nona[,3]))
genus_repro_hormones_corrected = bind_cols(genus_repro_hormones_nona, 
                                           as.data.frame(fdrtool(genus_repro_hormones_nona[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(genus_repro_hormones_corrected, "genus_repro_hormones_pa_corrected.csv")

genus_rainfall_hormones = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                                    as.data.frame(rainfall_pvaluematrix[,2]))
write.csv(genus_rainfall_hormones, "genus_rainfall_hormones_pa.csv")
genus_rainfall_hormones_nona = filter(genus_rainfall_hormones, 
                                      genus_rainfall_hormones[,3] != 0)
genus_rainfall_hormones_nona[,3] = as.numeric(as.character(genus_rainfall_hormones_nona[,3]))
genus_rainfall_hormones_corrected = bind_cols(genus_rainfall_hormones_nona, 
                                              as.data.frame(fdrtool(genus_rainfall_hormones_nona[,3], 
                                                                    statistic = "pvalue", plot = F)))
write.csv(genus_rainfall_hormones_corrected, "genus_rainfall_hormones_pa_corrected.csv")

genus_prog_hormones = bind_cols(as.data.frame(p_estimatematrix[,1:2]), 
                                as.data.frame(p_pvaluematrix[,2]))
write.csv(genus_prog_hormones, "genus_prog_hormones_pa.csv")
genus_prog_hormones_nona = filter(genus_prog_hormones, 
                                  genus_prog_hormones[,3] != 0)
genus_prog_hormones_nona[,3] = as.numeric(as.character(genus_prog_hormones_nona[,3]))
genus_prog_hormones_corrected = bind_cols(genus_prog_hormones_nona, 
                                          as.data.frame(fdrtool(genus_prog_hormones_nona[,3], 
                                                                statistic = "pvalue", plot = F)))
write.csv(genus_prog_hormones_corrected, "genus_prog_hormones_pa_corrected.csv")

genus_estr_hormones = bind_cols(as.data.frame(e_estimatematrix[,1:2]), 
                                as.data.frame(e_pvaluematrix[,2]))
write.csv(genus_estr_hormones, "genus_estr_hormones_pa.csv")
genus_estr_hormones_nona = filter(genus_estr_hormones, 
                                  genus_estr_hormones[,3] != 0)
genus_estr_hormones_nona[,3] = as.numeric(as.character(genus_estr_hormones_nona[,3]))
genus_estr_hormones_corrected = bind_cols(genus_estr_hormones_nona, 
                                          as.data.frame(fdrtool(genus_estr_hormones_nona[,3], 
                                                                statistic = "pvalue", plot = F)))
write.csv(genus_estr_hormones_corrected, "genus_estr_hormones_pa_corrected.csv")

genus_pperiod_hormones = bind_cols(as.data.frame(pp_estimatematrix[,1:2]), 
                                   as.data.frame(pp_pvaluematrix[,2]))
write.csv(genus_pperiod_hormones, "genus_pperiod_hormones_pa.csv")
genus_pperiod_hormones_nona = filter(genus_pperiod_hormones, 
                                     genus_pperiod_hormones[,3] != "NaN")
genus_pperiod_hormones_nona[,3] = as.numeric(as.character(genus_pperiod_hormones_nona[,3]))
genus_pperiod_hormones_corrected = bind_cols(genus_pperiod_hormones_nona, 
                                             as.data.frame(fdrtool(genus_pperiod_hormones_nona[,3], 
                                                                   statistic = "pvalue", plot = F)))
write.csv(genus_pperiod_hormones_corrected, "genus_pperiod_hormones_pa_corrected.csv")

genus_ebypp_hormones = bind_cols(as.data.frame(ebypp_estimatematrix[,1:2]), 
                                 as.data.frame(ebypp_pvaluematrix[,2]))
write.csv(genus_ebypp_hormones, "genus_ebypp_hormones_pa.csv")
genus_ebypp_hormones_nona = filter(genus_ebypp_hormones, 
                                   genus_ebypp_hormones[,3] != "NaN")
genus_ebypp_hormones_nona[,3] = as.numeric(as.character(genus_ebypp_hormones_nona[,3]))
genus_ebypp_hormones_corrected = bind_cols(genus_ebypp_hormones_nona, 
                                           as.data.frame(fdrtool(genus_ebypp_hormones_nona[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(genus_ebypp_hormones_corrected, "genus_ebypp_hormones_pa_corrected.csv")

#PERMANOVAs on unrarefied data----

#Import data, sort by sampleid prior to importing
unweighted_un = as.dist(read.table("unweighted-unrarefied-distance-matrix.tsv", header = T))
weighted_un = as.dist(read.table("weighted-unrarefied-distance-matrix.tsv", header = T))
unweighted_un_h = as.dist(read.table("unweighted-unrarefied-hormones-distance-matrix.tsv", header = T))
weighted_un_h = as.dist(read.table("weighted-unrarefied-hormones-distance-matrix.tsv", header = T))

metadata = read.csv("leaf_map_r_phyla.csv", header = T)
metadata_h = read.csv("leaf_map_r_allhormones.csv", header = T)

#PERMANOVAs
library(vegan)

adonis(unweighted_un~Repro_status + Rainfall, strata=metadata$Individual, 
       data=metadata, permutations=5000)
adonis(weighted_un~Repro_status + Rainfall, strata=metadata$Individual, 
       data=metadata, permutations=5000)
adonis(unweighted_un_h~logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
       strata=metadata_h$Individual, data=metadata_h, permutations=5000)
adonis(weighted_un_h~logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
       strata=metadata_h$Individual, data=metadata_h, permutations=5000)

##Normalization----

#Read in feature tables

library(tidyverse)

full = read.table("feature-table-full.tsv", sep = "\t", header = T)
level2 = read.table("feature-table-level2.tsv", sep = "\t", header = T)
level5 = read.table("feature-table-level5.tsv", sep = "\t", header = T)
level6 = read.table("feature-table-level6.tsv", sep = "\t", header = T)
level7 = read.table("feature-table-level7.tsv", sep = "\t", header = T)

meta = read.csv("leaf_map_r_phyla.csv", header = T)

#Create phyloseq objects

library(phyloseq)

metadata = sample_data(column_to_rownames(meta, "SampleID"))
full_otu= otu_table(column_to_rownames(full, "OTUID"), taxa_are_rows = T)
full_phylo = phyloseq(full_otu, metadata)
level2_otu = otu_table(column_to_rownames(level2, "OTUID"), taxa_are_rows = T)
level2_phylo = phyloseq(level2_otu, metadata)
level5_otu = otu_table(column_to_rownames(level5, "OTUID"), taxa_are_rows = T)
level5_phylo = phyloseq(level5_otu, metadata)
level6_otu = otu_table(column_to_rownames(level6, "OTUID"), taxa_are_rows = T)
level6_phylo = phyloseq(level6_otu, metadata)

#Normalize counts using variance stabilizing transformation

library(DESeq2)

full_deseq = phyloseq_to_deseq2(full_phylo, ~Repro_status)
full_deseq = estimateSizeFactors(full_deseq, type = "poscounts")
full_deseq_vst = varianceStabilizingTransformation(full_deseq, blind = FALSE, 
                                                   fitType = "parametric")
full_deseq_vst_out = as.data.frame(assay(full_deseq_vst)) %>% rownames_to_column()
deseq2_full = as.data.frame(t(column_to_rownames(full_deseq_vst_out))) %>% 
  rownames_to_column(var = "SampleID")
full = inner_join(deseq2_full, meta, by = "SampleID")
write.csv(full, "feature_full_deseq_vst.csv")

level2_deseq = phyloseq_to_deseq2(level2_phylo, ~Repro_status)
level2_deseq_vst = varianceStabilizingTransformation(level2_deseq, blind = FALSE, 
                                                     fitType = "parametric")
level2_deseq_vst_out = as.data.frame(assay(level2_deseq_vst)) %>% rownames_to_column()
deseq2_phyla = as.data.frame(t(column_to_rownames(level2_deseq_vst_out))) %>% 
  rownames_to_column(var = "SampleID")
phyla = inner_join(meta, deseq2_phyla, by = "SampleID")
write.csv(phyla, "feature_phyla_deseq_vst.csv")

level5_deseq = phyloseq_to_deseq2(level5_phylo, ~Repro_status)
level5_deseq_vst = varianceStabilizingTransformation(level5_deseq, blind = FALSE, 
                                                     fitType = "parametric")
level5_deseq_vst_out = as.data.frame(assay(level5_deseq_vst)) %>% rownames_to_column()
deseq2_family = as.data.frame(t(column_to_rownames(level5_deseq_vst_out))) %>% 
  rownames_to_column(var = "SampleID")
family = inner_join(meta, deseq2_family, by = "SampleID")
write.csv(family, "feature_family_deseq_vst.csv")

level6_deseq = phyloseq_to_deseq2(level6_phylo, ~Repro_status)
level6_deseq_vst = varianceStabilizingTransformation(level6_deseq, blind = FALSE, 
                                                     fitType = "parametric")
level6_deseq_vst_out = as.data.frame(assay(level6_deseq_vst)) %>% rownames_to_column()
deseq2_genus = as.data.frame(t(column_to_rownames(level6_deseq_vst_out))) %>% 
  rownames_to_column(var = "SampleID")
genus = inner_join(meta, deseq2_genus, by = "SampleID")
write.csv(genus, "feature_genus_deseq_vst.csv")

##LMEs for normalized count data----

library(nlme)
library(fdrtool)

#Phyla

acti = lme(fixed=Actinobacteria ~ Repro_status + Rainfall, 
           data = phyla, random = ~1|Individual)
anova(acti)

acti = lme(fixed=Actinobacteria ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
           data=phyla, random= ~1|Individual, na.action = na.omit)
anova(acti)

bact = lme(fixed=Bacteroidetes ~ Repro_status + Rainfall, 
           data = phyla, random = ~1|Individual)
anova(bact)

bact = lme(fixed=Bacteroidetes ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
           data=phyla, random= ~1|Individual, na.action = na.omit)
anova(bact)

firm = lme(fixed=Firmicutes ~ Repro_status + Rainfall, 
           data = phyla, random = ~1|Individual)
anova(firm)

firm = lme(fixed=Firmicutes ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
           data=phyla, random= ~1|Individual, na.action = na.omit)
anova(firm)

planc = lme(fixed=Planctomycetes ~ Repro_status + Rainfall, 
            data = phyla, random = ~1|Individual)
anova(planc)

planc = lme(fixed=Planctomycetes ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
            data=phyla, random= ~1|Individual, na.action = na.omit)
anova(planc)

prot = lme(fixed=Proteobacteria ~ Repro_status + Rainfall, 
           data = phyla, random = ~1|Individual)
anova(prot)

prot = lme(fixed=Proteobacteria ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
           data=phyla, random= ~1|Individual, na.action = na.omit)
anova(prot)

spiro = lme(fixed=Spirochaetes ~ Repro_status + Rainfall, 
            data = phyla, random = ~1|Individual)
anova(spiro)

spiro = lme(fixed=Spirochaetes ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
            data=phyla, random= ~1|Individual, na.action = na.omit)
anova(spiro)

tene = lme(fixed=Tenericutes ~ Repro_status + Rainfall, 
           data = phyla, random = ~1|Individual)
anova(tene)

tene = lme(fixed=Tenericutes ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
           data=phyla, random= ~1|Individual, na.action = na.omit)
anova(tene)

#Family
repro_estimatematrix = mat.or.vec(128,2)
repro_pvaluematrix = mat.or.vec(128,2)
rainfall_estimatematrix = mat.or.vec(128,2)
rainfall_pvaluematrix = mat.or.vec(128,2)

for(i in 21:147) {
  variable = family[,i]
  b<-lme(fixed=variable~Repro_status + Rainfall, data=family, 
         random= ~1|Individual)
  anova = anova(b)
  repro_estimatematrix[i-20,1] = names(family)[i]
  repro_estimatematrix[i-20,2] = anova[2,3]
  repro_pvaluematrix[i-20] = anova[2,4]
  rainfall_estimatematrix[i-20,1] = names(family)[i]
  rainfall_estimatematrix[i-20,2] = anova[3,3]
  rainfall_pvaluematrix[i-20] = anova[3,4]
}

family_repro = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                         as.data.frame(repro_pvaluematrix[,1]))
write.csv(family_repro, "family_repro_deseq.csv")
family_repro_corrected = bind_cols(family_repro_nona, 
                                   as.data.frame(fdrtool(family_repro[,3], 
                                                         statistic = "pvalue", plot = F)))
write.csv(family_repro_corrected, "family_repro_deseq_corrected.csv")

family_rainfall = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                            as.data.frame(rainfall_pvaluematrix[,1]))
write.csv(family_rainfall, "family_rainfall_deseq.csv")
family_rainfall_corrected = bind_cols(family_rainfall, 
                                      as.data.frame(fdrtool(family_rainfall_nona[,3], 
                                                            statistic = "pvalue", plot = F)))
write.csv(family_rainfall_corrected, "family_rainfall_deseq_corrected.csv")

family_hormones = family %>% filter(logP_avg3 > 0)

repro_estimatematrix = mat.or.vec(128,2)
repro_pvaluematrix = mat.or.vec(128,2)
rainfall_estimatematrix = mat.or.vec(128,2)
rainfall_pvaluematrix = mat.or.vec(128,2)
p_estimatematrix = mat.or.vec(128,2)
p_pvaluematrix = mat.or.vec(128,2)
e_estimatematrix = mat.or.vec(128,2)
e_pvaluematrix = mat.or.vec(128,2)
pperiod_estimatematrix = mat.or.vec(128,2)
pperiod_pvaluematrix = mat.or.vec(128,2)
pperiodbyfP_estimatematrix = mat.or.vec(128,2)
pperiodbyfP_pvaluematrix = mat.or.vec(128,2)

for(i in 21:147) {
  variable = family_hormones[,i]
  b = lme(fixed = variable ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
          data=family_hormones, random= ~1|Individual, control = lmeControl(returnObject = TRUE))
  anova = anova(b)
  repro_estimatematrix[i-20,1] = names(family_hormones)[i]
  repro_estimatematrix[i-20,2] = anova[5,3]
  repro_pvaluematrix[i-20] = anova[5,4]
  rainfall_estimatematrix[i-20,1] = names(family_hormones)[i]
  rainfall_estimatematrix[i-20,2] = anova[6,3]
  rainfall_pvaluematrix[i-20] = anova[6,4]
  p_estimatematrix[i-20,1] = names(family_hormones)[i]
  p_estimatematrix[i-20,2] = anova[2,3]
  p_pvaluematrix[i-20] = anova[2,4]
  e_estimatematrix[i-20,1] = names(family_hormones)[i]
  e_estimatematrix[i-20,2] = anova[4,3]
  e_pvaluematrix[i-20] = anova[4,4]
  pperiod_estimatematrix[i-20,1] = names(family_hormones)[i]
  pperiod_estimatematrix[i-20,2] = anova[3,3]
  pperiod_pvaluematrix[i-20] = anova[3,4]
  pperiodbyfP_estimatematrix[i-20,1] = names(family_hormones)[i]
  pperiodbyfP_estimatematrix[i-20,2] = anova[7,3]
  pperiodbyfP_pvaluematrix[i-20] = anova[7,4]
}

family_repro_hormones = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                         as.data.frame(repro_pvaluematrix[,1]))
write.csv(family_repro_hormones, "family_repro_hormones_deseq.csv")
family_repro_hormones_corrected = bind_cols(family_repro_hormones, 
                                   as.data.frame(fdrtool(family_repro_hormones[,3], 
                                                         statistic = "pvalue", plot = F)))
write.csv(family_repro_hormones_corrected, "family_repro_hormones_deseq_corrected.csv")

family_rainfall_hormones = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                                  as.data.frame(rainfall_pvaluematrix[,1]))
write.csv(family_rainfall_hormones, "family_rainfall_hormones_deseq.csv")
family_rainfall_hormones_corrected = bind_cols(family_rainfall_hormones, 
                                            as.data.frame(fdrtool(family_rainfall_hormones[,3], 
                                                                  statistic = "pvalue", plot = F)))
write.csv(family_rainfall_hormones_corrected, "family_rainfall_hormones_deseq_corrected.csv")

family_estr_hormones = bind_cols(as.data.frame(e_estimatematrix[,1:2]), 
                                     as.data.frame(e_pvaluematrix[,1]))
write.csv(family_estr_hormones, "family_estr_hormones_deseq.csv")
family_estr_hormones_corrected = bind_cols(family_estr_hormones, 
                                               as.data.frame(fdrtool(family_estr_hormones[,3], 
                                                                     statistic = "pvalue", plot = F)))
write.csv(family_estr_hormones_corrected, "family_estr_hormones_deseq_corrected.csv")

family_prog_hormones = bind_cols(as.data.frame(p_estimatematrix[,1:2]), 
                                 as.data.frame(p_pvaluematrix[,1]))
write.csv(family_prog_hormones, "family_prog_hormones_deseq.csv")
family_prog_hormones_corrected = bind_cols(family_prog_hormones, 
                                           as.data.frame(fdrtool(family_prog_hormones[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(family_prog_hormones_corrected, "family_prog_hormones_deseq_corrected.csv")

family_pperiod_hormones = bind_cols(as.data.frame(pperiod_estimatematrix[,1:2]), 
                                 as.data.frame(pperiod_pvaluematrix[,1]))
write.csv(family_pperiod_hormones, "family_pperiod_hormones_deseq.csv")
family_pperiod_hormones_corrected = bind_cols(family_pperiod_hormones, 
                                           as.data.frame(fdrtool(family_pperiod_hormones[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(family_pperiod_hormones_corrected, "family_pperiod_hormones_deseq_corrected.csv")

family_pperiodbyfP_hormones = bind_cols(as.data.frame(pperiodbyfP_estimatematrix[,1:2]), 
                                    as.data.frame(pperiod_pvaluematrix[,1]))
write.csv(family_pperiodbyfP_hormones, "family_pperiodbyfP_hormones_deseq.csv")
family_pperiodbyfP_hormones_corrected = bind_cols(family_pperiodbyfP_hormones, 
                                              as.data.frame(fdrtool(family_pperiodbyfP_hormones[,3], 
                                                                    statistic = "pvalue", plot = F)))
write.csv(family_pperiodbyfP_hormones_corrected, "family_pperiodbyfP_hormones_deseq_corrected.csv")

#Genus
repro_estimatematrix = mat.or.vec(185,2)
repro_pvaluematrix = mat.or.vec(185,2)
rainfall_estimatematrix = mat.or.vec(185,2)
rainfall_pvaluematrix = mat.or.vec(185,2)

for(i in 21:204) {
  variable = genus[,i]
  b<-lme(fixed=variable~Repro_status + Rainfall, data=genus, 
         random= ~1|Individual)
  anova = anova(b)
  repro_estimatematrix[i-20,1] = names(genus)[i]
  repro_estimatematrix[i-20,2] = anova[2,3]
  repro_pvaluematrix[i-20] = anova[2,4]
  rainfall_estimatematrix[i-20,1] = names(genus)[i]
  rainfall_estimatematrix[i-20,2] = anova[3,3]
  rainfall_pvaluematrix[i-20] = anova[3,4]
}

genus_repro = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                         as.data.frame(repro_pvaluematrix[,1]))
write.csv(genus_repro, "genus_repro_deseq.csv")
genus_repro_corrected = bind_cols(genus_repro, 
                                   as.data.frame(fdrtool(genus_repro[,3], 
                                                         statistic = "pvalue", plot = F)))
write.csv(genus_repro_corrected, "genus_repro_deseq_corrected.csv")

genus_rainfall = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                            as.data.frame(rainfall_pvaluematrix[,1]))
write.csv(genus_rainfall, "genus_rainfall_deseq.csv")
genus_rainfall_corrected = bind_cols(genus_rainfall, 
                                      as.data.frame(fdrtool(genus_rainfall[,3], 
                                                            statistic = "pvalue", plot = F)))
write.csv(genus_rainfall_corrected, "genus_rainfall_deseq_corrected.csv")

genus_hormones = genus %>% filter(logP_avg3 > 0)

repro_estimatematrix = mat.or.vec(185,2)
repro_pvaluematrix = mat.or.vec(185,2)
rainfall_estimatematrix = mat.or.vec(185,2)
rainfall_pvaluematrix = mat.or.vec(185,2)
p_estimatematrix = mat.or.vec(185,2)
p_pvaluematrix = mat.or.vec(185,2)
e_estimatematrix = mat.or.vec(185,2)
e_pvaluematrix = mat.or.vec(185,2)
pperiod_estimatematrix = mat.or.vec(185,2)
pperiod_pvaluematrix = mat.or.vec(185,2)
pperiodbyfP_estimatematrix = mat.or.vec(185,2)
pperiodbyfP_pvaluematrix = mat.or.vec(185,2)

for(i in 21:204) {
  variable = genus_hormones[,i]
  b = lme(fixed = variable ~ logP_avg3*Pperiod + logE_avg3 + Repro_status + Rainfall, 
          data=genus_hormones, random= ~1|Individual, control = lmeControl(returnObject = TRUE))
  anova = anova(b)
  repro_estimatematrix[i-20,1] = names(genus_hormones)[i]
  repro_estimatematrix[i-20,2] = anova[5,3]
  repro_pvaluematrix[i-20] = anova[5,4]
  rainfall_estimatematrix[i-20,1] = names(genus_hormones)[i]
  rainfall_estimatematrix[i-20,2] = anova[6,3]
  rainfall_pvaluematrix[i-20] = anova[6,4]
  p_estimatematrix[i-20,1] = names(genus_hormones)[i]
  p_estimatematrix[i-20,2] = anova[2,3]
  p_pvaluematrix[i-20] = anova[2,4]
  e_estimatematrix[i-20,1] = names(genus_hormones)[i]
  e_estimatematrix[i-20,2] = anova[4,3]
  e_pvaluematrix[i-20] = anova[4,4]
  pperiod_estimatematrix[i-20,1] = names(genus_hormones)[i]
  pperiod_estimatematrix[i-20,2] = anova[3,3]
  pperiod_pvaluematrix[i-20] = anova[3,4]
  pperiodbyfP_estimatematrix[i-20,1] = names(genus_hormones)[i]
  pperiodbyfP_estimatematrix[i-20,2] = anova[7,3]
  pperiodbyfP_pvaluematrix[i-20] = anova[7,4]
}

genus_repro_hormones = bind_cols(as.data.frame(repro_estimatematrix[,1:2]), 
                                  as.data.frame(repro_pvaluematrix[,1]))
write.csv(genus_repro_hormones, "genus_repro_hormones_deseq.csv")
genus_repro_hormones_corrected = bind_cols(genus_repro_hormones, 
                                            as.data.frame(fdrtool(genus_repro_hormones[,3], 
                                                                  statistic = "pvalue", plot = F)))
write.csv(genus_repro_hormones_corrected, "genus_repro_hormones_deseq_corrected.csv")

genus_rainfall_hormones = bind_cols(as.data.frame(rainfall_estimatematrix[,1:2]), 
                                     as.data.frame(rainfall_pvaluematrix[,1]))
write.csv(genus_rainfall_hormones, "genus_rainfall_hormones_deseq.csv")
genus_rainfall_hormones_corrected = bind_cols(genus_rainfall_hormones, 
                                               as.data.frame(fdrtool(genus_rainfall_hormones[,3], 
                                                                     statistic = "pvalue", plot = F)))
write.csv(genus_rainfall_hormones_corrected, "genus_rainfall_hormones_deseq_corrected.csv")

genus_estr_hormones = bind_cols(as.data.frame(e_estimatematrix[,1:2]), 
                                 as.data.frame(e_pvaluematrix[,1]))
write.csv(genus_estr_hormones, "genus_estr_hormones_deseq.csv")
genus_estr_hormones_corrected = bind_cols(genus_estr_hormones, 
                                           as.data.frame(fdrtool(genus_estr_hormones[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(genus_estr_hormones_corrected, "genus_estr_hormones_deseq_corrected.csv")

genus_prog_hormones = bind_cols(as.data.frame(p_estimatematrix[,1:2]), 
                                 as.data.frame(p_pvaluematrix[,1]))
write.csv(genus_prog_hormones, "genus_prog_hormones_deseq.csv")
genus_prog_hormones_corrected = bind_cols(genus_prog_hormones, 
                                           as.data.frame(fdrtool(genus_prog_hormones[,3], 
                                                                 statistic = "pvalue", plot = F)))
write.csv(genus_prog_hormones_corrected, "genus_prog_hormones_deseq_corrected.csv")

genus_pperiod_hormones = bind_cols(as.data.frame(pperiod_estimatematrix[,1:2]), 
                                    as.data.frame(pperiod_pvaluematrix[,1]))
write.csv(genus_pperiod_hormones, "genus_pperiod_hormones_deseq.csv")
genus_pperiod_hormones_corrected = bind_cols(genus_pperiod_hormones, 
                                              as.data.frame(fdrtool(genus_pperiod_hormones[,3], 
                                                                    statistic = "pvalue", plot = F)))
write.csv(genus_pperiod_hormones_corrected, "genus_pperiod_hormones_deseq_corrected.csv")

genus_pperiodbyfP_hormones = bind_cols(as.data.frame(pperiodbyfP_estimatematrix[,1:2]), 
                                        as.data.frame(pperiod_pvaluematrix[,1]))
write.csv(genus_pperiodbyfP_hormones, "genus_pperiodbyfP_hormones_deseq.csv")
genus_pperiodbyfP_hormones_corrected = bind_cols(genus_pperiodbyfP_hormones, 
                                                  as.data.frame(fdrtool(genus_pperiodbyfP_hormones[,3], 
                                                                        statistic = "pvalue", plot = F)))
write.csv(genus_pperiodbyfP_hormones_corrected, "genus_pperiodbyfP_hormones_deseq_corrected.csv")




