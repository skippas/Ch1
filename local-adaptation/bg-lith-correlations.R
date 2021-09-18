library(tidyverse)
library(ggpmisc)
library(ggrepel)
source("loading-cleaning.R")

# Theme ####
theme_set(theme_bw())
theme_update(legend.position = 'none')

# Plotting df ####
popMeans <- df %>%
  group_by(collapseSpp,substrate, abbrevs) %>%
  summarise(lumMean = mean(lumMean),
            xCoordMean = mean(xCoord),
            yCoordMean = mean(yCoord)) %>%
  pivot_wider(names_from = substrate,
              values_from = c(lumMean,xCoordMean, yCoordMean)) %>% ungroup() %>%
  mutate(siteLabs = unlist(lapply(str_split(abbrevs, "_"), `[[`, 1)))

# correlation lith rock ?means? ####
#* Overall ####
theme_set(theme_bw())
theme_update(aspect.ratio = 1,
            legend.position = 'none',
            plot.margin = unit(c(0,0,0,0), 'cm'))

pLumOvrl <- ggplot(popMeans, aes(x = lumMean_b, y = lumMean_a))+
  geom_point(aes(color = collapseSpp))+
  geom_abline(intercept = 0, lty = 2)+
  coord_fixed(xlim = c(.05,.5), ylim = c(.05,.5))+
  geom_smooth(method = "lm", color = "black", formula = y~x)+
  xlab("Rock luminance channel mean")+
  ylab("Lithops luminance channel mean")+
  stat_poly_eq(formula = y~x, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                                sep = "~~~")), parse = TRUE) 
png("output//colour-correlation//LR-lum.png", type = 'cairo',
    width = 215, height = 140, units = "mm", res = 300, pointsize = 6) 
pLumOvrl; dev.off()

pGROvrl <- ggplot(popMeans, aes(x = xCoordMean_b, y = xCoordMean_a))+
  geom_point(aes(color = collapseSpp))+
  geom_abline(intercept = 0, lty = 2)+
  coord_fixed(xlim = c(), ylim = c())+
  geom_smooth(method = "lm", color = "black", formula = y~x)+
  xlab("Rock G:R opponent channel mean")+ 
  ylab("Lithops G:R opponent channel mean")+
  stat_poly_eq(formula = y~x, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                                sep = "~~~")), parse = TRUE) 
png("output//colour-correlation//LR-GreenRed-correlation.png", type = 'cairo',
    width = 215,height = 140, units = "mm", res = 300, pointsize = 6)
pGROvrl; dev.off()

pBYOvrl <- ggplot(popMeans, aes(x = yCoordMean_b, y = yCoordMean_a))+
  geom_point(aes(color = collapseSpp))+
  geom_abline(intercept = 0, lty = 2)+
  coord_fixed(xlim = c(-1.6, -.35), ylim = c(-1.6, -.45))+
  geom_smooth(method = "lm", color = "black", formula = y~x)+
  xlab("Rock B:Y opponent channel mean")+
  ylab("Lithops B:Y opponent channel mean")+
  stat_poly_eq(formula = y~x, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                                sep = "~~~")), parse = TRUE)
png("output//colour-correlation//LR-BlueYell-correlation.png", type = 'cairo',
    width = 215, height = 140, units = "mm", res = 300, pointsize = 6) 
pBYOvrl; dev.off()


#* facet by species ####

pLumFacets <- ggplot(popMeans, aes(x = lumMean_b, y = lumMean_a)) +
  geom_point(data = popMeans[, c("lumMean_a", "lumMean_b")], color = "gray")+
  geom_point(aes(fill = collapseSpp), size = 2, pch =21)+ 
  geom_abline(intercept = 0, lty = 2)+
  geom_text_repel(aes(label = siteLabs))+
#  xlab("Rock luminance channel mean")+
#  ylab("Lithops luminance channel mean")+
  xlab("")+  ylab("")+
  facet_wrap(vars(collapseSpp), ncol = 5)+
  theme(legend.position = 'none', axis.title = element_blank())
png("output//colour-correlation//LR-lum-facetSpp-corr.png", type = 'cairo',
    width = 215, height = 140, units = "mm", res = 300, pointsize = 6)
pLumFacets; dev.off()

pGRFacets <- ggplot(popMeans, aes(x = xCoordMean_b, y = xCoordMean_a)) +
  geom_point(data = popMeans[, c("xCoordMean_a", "xCoordMean_b")], color = "gray")+
  geom_point(aes(fill = collapseSpp), size = 2, pch =21)+ 
  geom_abline(intercept = 0, lty = 2)+
  geom_text_repel(aes(label = siteLabs))+
#  xlab("Rock G:R opponent channel mean")+
#  ylab("Lithops G:R opponent channel mean")+
  coord_fixed(xlim = c(), ylim = c())+
  facet_wrap(vars(collapseSpp), ncol = 5)+
  theme(legend.position = 'none', axis.title = element_blank())
png("output//colour-correlation//LR-GreenReds-facetSpp-corr.png", type = 'cairo',
    width = 215, height = 140, units = "mm", res = 300, pointsize = 6)
pGRFacets; dev.off()

pBYFacets <- ggplot(popMeans, aes(x = yCoordMean_b, y = yCoordMean_a))+
  geom_point(data = popMeans[, c("yCoordMean_a", "yCoordMean_b")], color = "gray")+
  geom_point(aes(fill = collapseSpp), size = 2, pch =21)+ 
  geom_abline(intercept = 0, lty = 2)+
  geom_text_repel(aes(label = siteLabs))+
#  xlab("Rock B:Y opponent channel mean")+
#  ylab("Lithops B:Y opponent channel mean")+
  coord_fixed(xlim = c(), ylim = c())+
  facet_wrap(vars(collapseSpp), ncol = 5)+
  theme(legend.position = 'none', axis.title = element_blank())
png("output//colour-correlation//LR-blueYell-facetSpp-corr.png", type = 'cairo',
    width = 215, height = 140, units = "mm", res = 300, pointsize = 6) 
pBYFacets; dev.off()

#* Composite overall and species facets ####
library(patchwork)

png("output//colour-correlation//LR-lum-composite-corr.png", type = 'cairo',
    width = 215, height = 300, units = "mm", res = 300, pointsize = 6) 
pLumOvrl + pLumFacets + plot_layout(nrow = 2)
dev.off()

png("output//colour-correlation//LR-GR-composite-corr.png", type = 'cairo',
    width = 215, height = 300, units = "mm", res = 300, pointsize = 6) 
pGROvrl + pGRFacets + plot_layout(nrow = 2)
dev.off()

png("output//colour-correlation//LR-BY-composite-corr.png", type = 'cairo',
    width = 215, height = 300, units = "mm", res = 300, pointsize = 6) 
pBYOvrl + pBYFacets + plot_layout(nrow = 2)
dev.off()


# correlations for nearest 10% ####
# Generate distances ####
source("functions\\nearSubset.R") 
source("functions\\coldist_effic.R")
source("functions\\distToOrigDF.R")
source("functions\\gMeanDists.R")

df <- df %>% filter(substrate == "a") %>% group_by(abbrevs, mspec) %>% 
  sample_n(size = 1) %>% rbind(., df[df$substrate != "a",])

###* All dists ####
jnds_c <- colDistEffic(bg_subs = "c") %>% 
  filter(abbrevs.x == abbrevs.y)
jnds_b <- colDistEffic(bg_subs = "b") %>%
  filter(abbrevs.x == abbrevs.y)

###* Nearest dists ####

# Nearest subset plots ####
data <- list(jnds_b, jnds_b, jnds_c, jnds_c, jnds_b, jnds_b, jnds_c, jnds_c)
locOn <- rep(T, 8)
acrIm <- rep(F,8) 
vf <- c(rep("dL",4), rep("dS", 4))
qt <- rep(c(10,2), 4)
distSubsets <- Map(nearSubset,data,locOn,acrIm, vf, qt)
names(distSubsets) <- c("b_10%_dL","b_50%_dL", "c_10%_dL", "c_50%_dL",
                        "b_10%_dS", "b_50%_dS", "c_10%_dS", "c_50%_dS")

nrst10 <- distSubsets[-c(grep('50', names(distSubsets)))]
nrst10 <- lapply(nrst10, distToOrigDF)

rm(locOn,acrIm, vf, qt, data, distSubsets)

# Regression - Lithops coordinates predicted by substrate coordinates ####
# Problem - lumMean isn't the perceptually transformed version, nor is it log scaled etc.

###* plotting dfs ####
nrst10Tf <- nrst10 %>% 
  rbindlist(idcol = T) %>%
  merge(., df[, c(colnames(.)[colnames(.) %in% colnames(df)],
                  "xCoord", "yCoord")]) %>%
  group_by(.id,substrate,abbrevs, mspec) %>% 
  summarise_at(c("xCoord", "yCoord", "lumMean"), mean) %>%
  group_by(.id,substrate, abbrevs) %>% # this step not all that NB, highly correlated w straight mean
  summarise_at(c("xCoord", "yCoord", "lumMean"), mean) %>%
  separate(., col = .id,
           into = c("bgSubstrate","nrstSet", "visInf")) %>%
  pivot_longer(cols = c(xCoord, yCoord, lumMean),
               names_to = "visInfo") %>% 
  filter((!(visInf == "dL" & visInfo != "lumMean")) &
           !(visInf == "dS" & visInfo == "lumMean")) %>%
  select(-c("visInf")) %>%
  pivot_wider(names_from = visInfo,
              values_from = c(value)) %>% 
  group_by(nrstSet, substrate, abbrevs) %>% # To get rid of the extra a per pop
  summarise_at(c("xCoord", "yCoord", "lumMean"), mean) %>%
  mutate(substrate = as_factor(substrate),
         substrate = fct_recode(substrate,
                                Lithops = "a",
                                Rock = "b",
                                Soil = "c")) %>%
  pivot_wider(names_from = substrate,
              values_from = c(lumMean, xCoord, yCoord)) %>%
  merge(., unique(df[, c('abbrevs', 'collapseSpp')])) %>%
  mutate(siteLabs = unlist(lapply(str_split(abbrevs, "_"), `[[`, 1)))

# correlation nearest 10% plots ####

pLumOvrl <- ggplot(nrst10Tf, aes(x = lumMean_Rock, y = lumMean_Lithops))+
  geom_point(aes(color = collapseSpp))+
  geom_abline(intercept = 0, lty = 2)+
  geom_smooth(method = "lm", color = "black", formula = y~x)+
  xlab("Rock luminance channel mean")+
  ylab("Lithops luminance channel mean")+
  stat_poly_eq(formula = y~x, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                                sep = "~~~")), parse = TRUE) 

pLumFacets <- ggplot(nrst10Tf, aes(x = lumMean_Rock, y = lumMean_Lithops))+
  geom_point(data = nrst10Tf[, c("lumMean_Lithops", "lumMean_Rock")], color = "gray")+
  geom_point(aes(fill = collapseSpp), size = 2, pch =21)+ 
  geom_abline(intercept = 0, lty = 2)+
  geom_text_repel(aes(label = siteLabs))+
  #  xlab("Rock luminance channel mean")+
  #  ylab("Lithops luminance channel mean")+
  xlab("")+  ylab("")+
  facet_wrap(vars(collapseSpp), ncol = 5)+
  theme(legend.position = 'none', axis.title = element_blank())

png("output//colour-correlation//LR-lum-comp-10-corr.png", type = 'cairo',
    width = 215, height = 300, units = "mm", res = 300, pointsize = 6)
pLumOvrl + pLumFacets + plot_layout(nrow = 2) ; dev.off()

pGROvrl <- ggplot(nrst10Tf, aes(x = xCoord_Rock, y = xCoord_Lithops))+
  geom_point(aes(color = collapseSpp))+
  geom_abline(intercept = 0, lty = 2)+
  geom_smooth(method = "lm", color = "black", formula = y~x)+
  xlab("Rock G:R opponent mean")+
  ylab("Lithops G:R opponent mean")+
  stat_poly_eq(formula = y~x, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                                sep = "~~~")), parse = TRUE) 

pGRFacets <- ggplot(nrst10Tf, aes(x = xCoord_Rock, y = xCoord_Lithops))+
  geom_point(data = nrst10Tf[, c("xCoord_Lithops", "xCoord_Rock")],
             color = "gray")+
  geom_point(aes(fill = collapseSpp), size = 2, pch =21)+ 
  geom_abline(intercept = 0, lty = 2)+
  geom_text_repel(aes(label = siteLabs))+
  #  xlab("Rock luminance channel mean")+
  #  ylab("Lithops luminance channel mean")+
  xlab("")+  ylab("")+
  facet_wrap(vars(collapseSpp), ncol = 5)+
  theme(legend.position = 'none', axis.title = element_blank())
png("output//colour-correlation//LR-GR-comp-10-corr.png", type = 'cairo',
    width = 215, height = 300, units = "mm", res = 300, pointsize = 6)
pGROvrl + pGRFacets + plot_layout(nrow = 2) ; dev.off()

pBYOvrl <- ggplot(nrst10Tf, aes(x = yCoord_Rock, y = yCoord_Lithops))+
  geom_point(aes(color = collapseSpp))+
  geom_abline(intercept = 0, lty = 2)+
  geom_smooth(method = "lm", color = "black", formula = y~x)+
  xlab("Rock B:Y opponent mean")+
  ylab("Lithops B:Y opponent mean")+
  stat_poly_eq(formula = y~x, aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                                sep = "~~~")), parse = TRUE) 

pBYFacets <- ggplot(nrst10Tf, aes(x = yCoord_Rock, y = yCoord_Lithops))+
  geom_point(data = nrst10Tf[, c("yCoord_Lithops", "yCoord_Rock")],
             color = "gray")+
  geom_point(aes(fill = collapseSpp), size = 2, pch =21)+ 
  geom_abline(intercept = 0, lty = 2)+
  geom_text_repel(aes(label = siteLabs))+
  #  xlab("Rock luminance channel mean")+
  #  ylab("Lithops luminance channel mean")+
  xlab("")+  ylab("")+
  facet_wrap(vars(collapseSpp), ncol = 5)+
  theme(legend.position = 'none', axis.title = element_blank())
png("output//colour-correlation//LR-BY-comp-10-corr.png", type = 'cairo',
    width = 215, height = 300, units = "mm", res = 300, pointsize = 6)
pBYOvrl + pBYFacets + plot_layout(nrow = 2) ; dev.off()




# Model with species as predictor ####
modelDf <- nrst10Tf %>% 
  filter(collapseSpp != 'singlePops') %>%
  mutate(collapseSpp = as.factor(collapseSpp))

sppModel <- with(modelDf,lm(lumMean_Lithops ~  collapseSpp))
summary(sppModel)

