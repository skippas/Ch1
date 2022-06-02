# SETUP ####
library(ggpubr)
library(rstatix)
source("loading-cleaning.R")
source("functions\\nearSubset.R") 
source("functions\\coldist_effic.R")
source("functions\\distToOrigDF.R")

theme_set(theme_bw())
theme_update(panel.spacing = unit(1, "lines"),
             strip.placement = "outside",
             strip.background = element_rect(colour = "black")) # see riffomonas proj


# Generate distances ####
set.seed(42)
df <- df %>% filter(substrate == "a") %>% group_by(abbrevs, mspec) %>% 
  sample_n(size = 1) %>%
  rbind(., df[df$substrate != "a",])

###* All dists 
jnds_c <- colDistEffic(bg_subs = "c") %>% 
  filter(abbrevs.x == abbrevs.y)
jnds_b <- colDistEffic(bg_subs = "b") %>%
  filter(abbrevs.x == abbrevs.y)

###* Nearest dists
# nearest 10% each substrate OR 
# nearest n ROIs across all images for a single plant in each image
nrstSub_dL_c <- nearSubset(jnds_c, localOnly = T, acrossImage = F, 
                           visfeature = "dL", qtiles =  10) 
nrstSub_dL_b <- nearSubset(jnds_b, localOnly = T, acrossImage = F, 
                           visfeature = "dL", qtiles =  10) 
nrstSub_dS_c <- nearSubset(jnds_c, localOnly = T, acrossImage = F, 
                           visfeature = "dS", qtiles =  10) 
nrstSub_dS_b <- nearSubset(jnds_b, localOnly = T, acrossImage = F, 
                           visfeature = "dS", qtiles =  10) # IS ANYTHING DIFF HAPPENING WHEN I change T/F?

homogDF <- function(df,dropInfo,substrate){ # format where substrate and visInfo can be disentangled 
  visInfoType <- c("dL", "dS")              # and kept in same df
  keepInfo <- visInfoType[!visInfoType %in% dropInfo] # horribly verbose function? 
  df <- df[!names(df) %in% dropInfo]                  
  df$subs <- substrate
  df$visInfo <- keepInfo
  names(df)[names(df) == keepInfo] <- "distance"
  return(df)
}
nrstSub_dL_b <- homogDF(df = nrstSub_dL_b,dropInfo = "dS", substrate = "b")
nrstSub_dL_c <- homogDF(df = nrstSub_dL_c,dropInfo = "dS", substrate = "c")
nrstSub_dS_b <- homogDF(df = nrstSub_dS_b,dropInfo = "dL", substrate = "b")
nrstSub_dS_c <- homogDF(df = nrstSub_dS_c,dropInfo = "dL", substrate = "c")
nrstSub_bc <- rbind(nrstSub_dS_b, nrstSub_dS_c, nrstSub_dL_b, nrstSub_dL_c)

rm(nrstSub_dL_b, nrstSub_dL_c, nrstSub_dS_b, nrstSub_dS_c)
###* Paired distance metric ####

source('functions//gmean.R')

# Nearest

nrstSub_bc <- nrstSub_bc %>% filter(mspec.x == mspec.y) # this is what happens when set across image F anyway?

# extracting the nearest original data. actually random lithops, not nearest.

lumDists <- lapply(nrstLumLs, coldist,
                    n = c(1,16,32), weber = .05, weber.achro = .05, 
                    weber.ref = "longest", qcatch = "Qi", achromatic = TRUE)

chromDists <- lapply(nrstChromLs, coldist, 
                      n = c(1,16,32), weber = .05, weber.achro = .05, 
                      weber.ref = "longest", qcatch = "Qi", achromatic = TRUE)

lumDists <- rbindlist(lumDists, idcol= "names") %>%
  separate(., names, into = c("abbrevs", "mspec"), sep = "--") %>%
  select(-c(dS)) %>%
  filter(patch1 == 'a' | patch2 == 'a')

chromDists <- rbindlist(chromDists, idcol= "names") %>%
  separate(., names, into = c("abbrevs", "mspec"), sep = "--") %>%
  select(-c(dL)) %>%
  filter(patch1 == 'a' | patch2 == 'a')

perImageDists <- merge(chromDists, lumDists) %>%
  pivot_longer(cols = c('dS', 'dL'), 
               values_to = 'distance', names_to = 'visInfo') %>%
  rename(., substrate = 'patch2') %>%
  select(-c(patch1)) %>%
  pivot_wider(names_from = substrate, values_from = distance) %>%
  merge(., unique(df[,c('abbrevs', 'collapseSpp')])) %>%
  mutate(dataset = 'Nearest 10%')

rm(chromDists, lumDists, lstLumDists, nrstChromROIs,
   nrstChromLs, nrstLumLs, nrstLumROIs,
   chromDists, lumDists, nrstLithops)

# All
logQi <- df %>% 
  group_by(abbrevs, mspec, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  mutate(splitFct = paste(abbrevs,mspec, sep = '--')) %>%
  as.data.frame() %>% ungroup() %>%
  split(., .$splitFct)

logQi <- lapply(logQi, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

logQi <- Filter(function(x) nrow(x)>=2, logQi)

perImageDistsAll <- lapply(logQi, coldist,
                           n = c(1,16,32), weber = .05, weber.achro = .05, 
                           weber.ref = "longest", qcatch = "Qi", achromatic = TRUE)

perImageDistsAll <- rbindlist(perImageDistsAll, idcol= "names") %>%
  separate(., names, into = c("abbrevs", "mspec"), sep = "--") %>%
  filter(patch1 == 'a' | patch2 == 'a') %>%
  merge(., unique(df[,c('abbrevs', 'collapseSpp')]),
        by.x = 'abbrevs', by.y = 'abbrevs') %>%
  mutate(dataset = 'Average') %>%
  pivot_longer(cols = c(dS,dL), 
               names_to = 'visInfo', values_to = 'distance') %>%
  pivot_wider(names_from = patch2, values_from = distance) %>%
  select(-c(patch1))

perImageDists <- rbind(perImageDists, perImageDistsAll) %>%
  mutate(facetVar = paste(visInfo, dataset))

rm(logQi, perImageDistsAll)

###** Paired T-test, anno df, boxplots
t_test <- function(df, mu = 0, alt = "two.sided", paired = T,
                   conf.level = .95,var.equal = F){
  tidy(t.test(df$b, df$c,
              mu = mu, alt = alt,
              conf.level = conf.level,
              paired = paired, var.equal = var.equal))
}

pairedTests <- perImageDists %>%
  group_by(abbrevs, visInfo, dataset) %>%
  nest() %>%
  mutate(ttest = map(data, t_test)) %>%
  unnest(ttest) %>%
  group_by(visInfo, dataset) %>%
  mutate(p.adj = p.adjust(p.value, n = 56, method = 'holm')) %>%
  merge(., unique(df[, c('abbrevs', 'collapseSpp')])) %>%
        mutate(labelz = case_when(p.adj < 0.05 ~ "*"),
               labelzPValue = case_when(p.value < 0.05 ~ '*'),
               dPerLith =
                 case_when(visInfo == 'dL' & dataset == 'Average' ~ 45,
                           visInfo == 'dL' & dataset == 'Nearest 10%' ~40, 
                           visInfo == 'dS' & dataset == 'Average' ~ 7,
                           TRUE ~ 7),
               significanceDirection = case_when(estimate > 0 ~ 'Soil< Rock',
                                                 TRUE ~ 'Rock< Soil'),
               facetVar = paste(visInfo, dataset))

# Results 

pairedTests %>% 
  filter(p.adj < 0.05) %>%
  group_by(dataset, visInfo) %>%
  summarise(n = n(),
            'R<S' = sum(statistic < 0),
            'S<R' = sum(statistic > 0),
            rockSigPcent = `R<S`/ n*100,
            soilSigPcent = `S<R` / n*100)

# mean JNDs < 1
inImageMeans %>% ungroup %>%
  filter(dataset == 'Average', visInfo == 'dL', c<1) %>%
  count()

###* boxlots
# boxplots idea stolen from https://www.stat.auckland.ac.nz/~paul/RGraphics/examples-dotplot.R

facetVarLabeller <- c('dL Nearest 10%' = 'Lum Nearest 10%',
                      'dS Nearest 10%' = 'Chroma Nearest 10%',
                      'dL Average' = 'Luminance Average',
                      'dS Average' = 'Chroma Average')
sppLabeller<- c(olivacea = 'L. oliv.', localis = 'L. loc.',
                marmorataSpp = 'L. marmorata', divergensSpp = 'L. divergens',
                halliiSpp = 'L. hallii', leslieiSpp = 'L. lesliei',
                otzeniana = 'L. otz.', compt. = 'L. compt.',
                fulleri = 'L. fulleri', bromf. = 'L. bromf.',
                hookeriSpp = 'L. hookeri', meyeri = 'L. me.', 
                dint. = 'L. din.', verruc. = 'L. verru.', 
                singlePops = 'singlePops')

p <- perImageDists %>% 
  pivot_longer(cols = c(b,c), 
               names_to = "subs",
               values_to = "dPerLith") %>%
  ggplot(aes(x = abbrevs, y =  dPerLith))+
  geom_boxplot(aes(fill = subs), outlier.shape = NA)+
  geom_text(data = pairedTests, key_glyph = 'point', # No clue how I got the legend to show a star
            aes(label = labelz, colour = significanceDirection,
                x = abbrevs, y = dPerLith), fontface = 'bold', size = 5)+
  geom_text(data = pairedTests, key_glyph = 'point', # No clue how I got the legend to show a star
            aes(label = labelzPValue, colour = significanceDirection,
                x = abbrevs, y = dPerLith), alpha =.4, size = 5)+
  geom_hline(yintercept = 1, lty = 2, colour = 'grey70')+
  geom_point(size = 0, stroke = 0)+
  scale_fill_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                    labels = c("lithops-rock", "lithops-soil"))+
  scale_colour_manual(values = c("Rock< Soil" = "#00BFC4",
                                 "Soil< Rock" = "#F8766D"))+
  guides(colour = guide_legend(
    override.aes = list(size = 5, shape = c(utf8ToInt("*"), utf8ToInt("*")))))+
  labs(fill = "Substrate contrast",
       colour = 'Significance direction',
       x = "Populations by species",
       y = "Noise scaled distances (JNDs)")+
  scale_y_continuous(sec.axis = dup_axis(name = waiver()))+
  coord_flip()+ 
  facet_grid(collapseSpp~facetVar,
             switch = "y", scales = "free", space = "free_y",
             labeller = labeller(facetVar = facetVarLabeller,
                                 collapseSpp = sppLabeller))+
  theme(strip.placement = 'outside',
        strip.text.y = element_text(face = 'italic'),
        strip.background = element_rect(fill = 'grey91'),
        panel.spacing.y = unit(0,'lines'),
        legend.position = 'top')

png(file = "output//rock-v-soil//RvS-pairedTest2.png", type = 'cairo',
    width = 215, height = 280, units = "mm", res = 300, pointsize = 6)
p
dev.off()

rm(pairedTests, p, perImageDists)

# LR V LS pop-mean distances (scatterplots) ####

nrstSub_bc <- nrstSub_bc %>% filter(mspec.x == mspec.y)

nrstLithops <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                   "abbrevs", "mspec","roi", "substrate")]), 
                     unique(nrstSub_bc[, c("abbrevs.x", "mspec.x", "roi.x")]),
                     by.y = c("abbrevs.x","mspec.x","roi.x"),
                     by.x = c("abbrevs","mspec", "roi")) %>%
  distinct(.) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) 

gMeanLum <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                   "abbrevs", "mspec","roi", "substrate")]), 
                     unique(nrstSub_bc[
                       nrstSub_bc$visInfo == 'dL',
                       c("abbrevs.x", "mspec.x", "roi.y")]),
                     by.y = c("abbrevs.x","mspec.x","roi.y"),
                     by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(abbrevs, substrate, mspec) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  rbind(nrstLithops) %>%
  as.data.frame(.) %>%
  split(.,.$abbrevs)

gMeanChrom <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                     "abbrevs", "mspec","roi", "substrate")]), 
                       unique(nrstSub_bc[
                         nrstSub_bc$visInfo == 'dS',
                         c("abbrevs.x", "mspec.x", "roi.y")]),
                       by.y = c("abbrevs.x","mspec.x","roi.y"),
                       by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(abbrevs, substrate, mspec) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  rbind(nrstLithops) %>%
  as.data.frame(.) %>%
  split(.,.$abbrevs)

gMeanAll <- df %>%  
  group_by(abbrevs, substrate, mspec) %>% # each image contributes equally to gmean dist
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>% 
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  as.data.frame(.) %>%
  split(.,.$abbrevs)

gMeanLum <- lapply(gMeanLum, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

gMeanChrom <- lapply(gMeanChrom, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

gMeanAll <- lapply(gMeanAll, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

gMeanLum <- lapply(gMeanLum, coldist,
                   n = c(1,16,32), weber = .05, weber.achro = .05, 
                   weber.ref = "longest",
                   qcatch = "Qi", achromatic = TRUE) %>%
  rbindlist(., idcol = 'abbrevs') %>%
  pivot_longer(cols = c(dS,dL), 
               names_to = 'visInfo', values_to = 'distance') %>%
  filter(patch1 == 'a' | patch2 == 'a', visInfo == 'dL') %>%
  select(-c(patch2)) %>%
  rename(substrate = 'patch1') %>%
  mutate(dataset = 'Nearest 10%')

gMeanChrom <- lapply(gMeanChrom, coldist,
                     n = c(1,16,32), weber = .05, weber.achro = .05, 
                     weber.ref = "longest",
                     qcatch = "Qi", achromatic = TRUE) %>%
  rbindlist(., idcol = 'abbrevs') %>%
  pivot_longer(cols = c(dS,dL), 
               names_to = 'visInfo', values_to = 'distance') %>%
  filter(patch1 == 'a' | patch2 == 'a', visInfo == 'dS') %>%
  select(-c(patch2)) %>%
  rename(substrate = 'patch1') %>%
  mutate(dataset = 'Nearest 10%')
  
gMeanAll <- lapply(gMeanAll, coldist,
                   n = c(1,16,32), weber = .05, weber.achro = .05, 
                   weber.ref = "longest",
                   qcatch = "Qi", achromatic = TRUE) %>%
  rbindlist(., idcol = 'abbrevs') %>%
  filter(patch1 == 'a' | patch2 == 'a') %>%
  select(-c(patch1)) %>%
  rename(substrate = 'patch2') %>%
  mutate(dataset = 'Average') %>%
  pivot_longer(cols = c(dS,dL), 
               names_to = 'visInfo', values_to = 'distance')

gMeans <- rbind(gMeanAll, gMeanLum, gMeanChrom) %>% 
  pivot_wider(names_from = 'substrate', values_from = 'distance') 
rm(gMeanAll, gMeanLum, gMeanChrom, nrstLithops)

# # Add substrate frequency variable
# source("rock-v-soil//rock-estimate.R")
# props <- df %>% count(abbrevs, substrate)
# props <- props %>% pivot_wider(names_from = substrate, values_from = n)
# props$rock_cover <- (props$b / (props$b + props$c)*100)
# props$soil_cover <- (props$c / (props$b + props$c)*100)
# props <- arrange(props, as.factor(abbrevs))
# allGmeanDist <- merge(props, allGmeanDist, by.x = "abbrevs", by.y = "bgpop")
# source("loading-cleaning.R") # not ideal. could change name of rock estimate df

p1<- ggplot(gMeans %>% filter(visInfo == "dL"), aes(x =c, y = b))+
  geom_point()+
  xlim(c(0,14))+ ylim(c(0,14))+
  xlab(NULL)+ylab('LR mean distance')+
  labs(title = "Luminance")+
  geom_abline(lty = 2)+
  geom_hline(yintercept = 1, lty = 2, colour = 'grey70')+
  geom_vline(xintercept = 1, lty = 2, colour = 'grey70')+
  coord_equal()+
  facet_grid(vars(dataset), switch = 'y')+
  theme(strip.placement = 'outside',
        strip.background = element_blank(),
        plot.title = element_text(hjust = .5, vjust = .5, size = 10),
        axis.title.x = element_text(
          margin = unit(c(0, 0, 0, 0), "mm")))

p2<- ggplot(gMeans %>% filter(visInfo == "dS"), aes(x =c, y = b))+
  geom_point()+
  xlim(c(0,6))+ ylim(c(0,6))+
  xlab(NULL)+ylab(NULL)+
  labs(title = "Chroma")+
  geom_abline(lty = 2)+
  geom_hline(yintercept = 1, lty = 2, colour = 'grey70')+
  geom_vline(xintercept = 1, lty = 2, colour = 'grey70')+
  coord_equal()+
  facet_grid(vars(dataset), switch = 'y')+
  theme(plot.title = element_text(hjust = .5, vjust = .5, size = 10),
        strip.text = element_blank(),
        axis.title.x = element_text(
          margin = unit(c(0, 0, 0, 0), "mm")))

library(patchwork)
png("output\\rock-v-soil\\LR-LS-distance-relship-varyNrst.png",
    type = 'cairo',units = 'mm', res = 300, width = 200, height = 200)
#cowplot::plot_grid(p1,p2, ncol = 2, align = 'v')
p1+p2+ plot_annotation(
  caption = 'LS mean distance',
  theme = theme(
    plot.caption = 
      element_text(hjust = .55,margin = margin(t=0), size = 11.5)))+
  plot_layout(ncol = 2)
dev.off()

rm(p,p1,p2)

# Results
# populations within 1 > JND of one or either substrate
gMeans %>%
  pivot_longer(cols = c(b,c),
               names_to = 'substrate',
               values_to = 'distance') %>%
  filter(distance < 1) %>%
  count(dataset, visInfo, substrate) %>%
  mutate(percentBelowSub = n/56*100)

gMeans %>%
  pivot_longer(cols = c(b,c),
               names_to = 'substrate',
               values_to = 'distance') %>%
  filter(distance < 1) %>%
  group_by(dataset, visInfo) %>%
  summarise(popsBelowEither = n_distinct(abbrevs)) %>%
  mutate(pcentBelowEither = popsBelowEither/56*100)

gMeans %>% 
  group_by(dataset, visInfo) %>%
  summarise(`R<S` = sum(b<c),
            `S<R` = sum(c<b),
            `S<R%` = `S<R`/56*100,
            `R<S%` = `R<S`/56*100)

###** Overall Boxplots pop-mean LR, LS dists ####
### paired t-test difference in overall LR-LS dists

overallTtest <- gMeans %>%
  group_by(dataset, visInfo) %>%
  nest() %>%
  mutate(ttest = map(data, t_test)) %>%
  unnest(ttest) %>% ungroup()

overallTtest <- overallTtest %>%
  mutate(labelz = 
           case_when(p.value < 0.001 ~ '***',
                     p.value < 0.01 ~ '**',
                     p.value < 0.05 ~ '*'))

pData <- gMeans %>%
  pivot_longer(cols = c("b","c"),
               names_to = "substrate", values_to = "meanDist")
p2 <-ggplot(pData, aes(dataset,meanDist))+
  geom_boxplot(aes(fill = substrate))+
  geom_text(data = overallTtest,
            aes(y = Inf,label = labelz), vjust = 1,size = 5)+
  theme_classic()+
  labs(y = 'Population-mean distances',
       fill = 'Substrate contrast')+
  scale_x_discrete('Dataset')+
  scale_fill_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                    labels = c("lithops-rock", "lithops-soil"))+
  facet_grid(vars(visInfo), scales = "free_y",
             labeller = labeller(visInfo = c('dL' = 'Luminance',
                                                'dS' = 'Chroma')))

png("output\\rock-v-soil\\LR-LS-dist-relship-varyNrstPcent-boxplots.png",
    type = 'cairo',units = 'mm', res = 300, width = 200, height = 140)
p2
dev.off()

# Results
rm(pData, overallTtest)

###** STRANGE RESULT ####
# TK verruc HG smaller distance L-S for the all subset because it has a 
# different subset of Lithops. Because there are only a few images with soils, 
# the nearest subsets calculate the location of only the few Lithops which are
# in those images  

###** Debugging ####

# # nearest subsets
# c_10_dL <- nearSubset(jnds_c, localOnly = T, acrossImage = F, 
#                            visfeature = "dL", qtiles =  10) 
# b_10_dL <- nearSubset(jnds_b, localOnly = T, acrossImage = F, 
#                            visfeature = "dL", qtiles =  10) 
# c_50_dL <- nearSubset(jnds_c, localOnly = T, acrossImage = F, 
#                       visfeature = "dL", qtiles =  2) 
# b_50_dL <- nearSubset(jnds_b, localOnly = T, acrossImage = F, 
#                       visfeature = "dL", qtiles =  2) 
# bc_50_dL <- rbind(c_50_dL, b_50_dL)
# bc_10_dL <- rbind(c_10_dL, b_10_dL)
# 
# lsNrst <- list(bc_50_dL, bc_10_dL)
# names(lsNrst) <- c("bc_50_dL", "bc_10_dL")
# 
# rm(c_10_dL,b_10_dL, b_50_dL, c_50_dL)
# 
# # coordinates of the nearest subsets
# lsNrst <- lapply(lsNrst, distToOrigDF)
# df2 <- df[,c(colnames(lsNrst[[1]]))]
# df2 <- ungroup(df2)
# lsNrst[["bc_all_dL"]] <- df2 # add in the 'all' dists
# 
# rm(df2)
# 
# db_distSubsets <- rbindlist(lsNrst, idcol = T)
# db_distSubsets <- db_distSubsets %>% 
#   group_by(.id, abbrevs, substrate, mspec) %>% 
#   summarise_at(vars(swMean, mwMean, lwMean, lumMean), mean) # each image contributes equally to gmean dist
# lsNrst <- split(db_distSubsets, db_distSubsets$.id)
# 
# # geometric mean distances between the coordinates
# lsNrst <- lapply(lsNrst, gMeanDists)
# db_gMeanSubsets <- rbindlist(lsNrst, idcol = T)
# db_gMeanSubsets <- db_gMeanSubsets %>% filter(comparison == "local", patch1 == "a") 
# db_gMeanSubsets <- db_gMeanSubsets %>% select(-c("dS")) %>% rename(gmDist = "dL")  # drop rows where visInf and visf !=
# 
# db_gMeanSubsets <- db_gMeanSubsets %>% select(-c("patch1")) %>% 
#   pivot_wider(values_from = "gmDist", names_from = "patch2")
# 
# ggplot(gMeanSubsets, aes(x =b, y = c))+
#   geom_point()+
#   #coord_equal()+
#   theme(aspect.ratio = 1)+
#   xlim(c(0,12))+ ylim(c(0,12))+
#   facet_grid(vars(.id), scales = 'free')
# 
# 
# rm(x,y,xx,yy)


# Regression - Lithops coordinates predicted by substrate coordinates ####

###* plotting dfs

nrstLithops <- merge(unique(df[, c('xCoord', 'yCoord',"lumCoord",
                                   "abbrevs", "mspec","roi", "substrate")]), 
                     unique(nrstSub_bc[, c("abbrevs.x", "mspec.x", "roi.x")]),
                     by.y = c("abbrevs.x","mspec.x","roi.x"),
                     by.x = c("abbrevs","mspec", "roi")) %>%
  distinct(.) %>%
  group_by(substrate,abbrevs) %>%
  summarise(xCoord = mean(xCoord),
            yCoord = mean(yCoord),
            lumCoord = mean(lumCoord)) %>%
  select(c(abbrevs,substrate, xCoord, yCoord, lumCoord))

gMeanLum <- merge(unique(df[, c("lumCoord",
                                "abbrevs", "mspec","roi", "substrate")]), 
                  unique(nrstSub_bc[
                    nrstSub_bc$visInfo == 'dL',
                    c("abbrevs.x", "mspec.x", "roi.y")]),
                  by.y = c("abbrevs.x","mspec.x","roi.y"),
                  by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(substrate,abbrevs, mspec) %>% 
  summarise(lumCoord = mean(lumCoord)) %>%
  group_by(substrate, abbrevs) %>% # this step not all that NB, highly correlated w straight mean
  summarise(lumCoord = mean(lumCoord)) 

meanCoordsChrom <- merge(unique(df[, c("xCoord", "yCoord",
                                  "abbrevs", "mspec","roi", "substrate")]), 
                    unique(nrstSub_bc[
                      nrstSub_bc$visInfo == 'dS',
                      c("abbrevs.x", "mspec.x", "roi.y")]),
                    by.y = c("abbrevs.x","mspec.x","roi.y"),
                    by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(abbrevs, substrate, mspec) %>% 
  summarise_at(c("xCoord", "yCoord"), mean) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(c("xCoord", "yCoord"), mean) %>%
  merge(., gMeanLum) %>%
  rbind(nrstLithops) %>%
  pivot_wider(names_from = substrate,
              values_from = c(xCoord, yCoord, lumCoord)) %>%
  pivot_longer(cols = c(contains(c("_b","_c"))), 
               names_to = c("predVisInfo", "predSub"),
               names_sep = "_",
               values_to = "predValue") %>% 
  pivot_longer(cols = c(contains(c("_a"))),
               names_to = c("lithVisInfo", "lithSub"),
               names_sep = "_",
               values_to = 'lithValue') %>%
  filter(predVisInfo == lithVisInfo) %>%
  select(-c(lithVisInfo, lithSub)) %>%
  mutate(dataset = "Nearest 10%")

meanAll <- df %>%  
  group_by(abbrevs, substrate, mspec) %>% 
  summarise(xCoord = mean(xCoord),
            yCoord = mean(yCoord),
            lumCoord = mean(lumCoord)) %>% 
  group_by(abbrevs, substrate) %>%
  summarise(xCoord = mean(xCoord),
            yCoord = mean(yCoord),
            lumCoord = mean(lumCoord)) %>% 
  pivot_wider(names_from = substrate,
              values_from = c(xCoord, yCoord, lumCoord)) %>%
  pivot_longer(cols = c(contains(c("_b","_c"))), 
               names_to = c("predVisInfo", "predSub"),
               names_sep = "_",
               values_to = "predValue") %>% 
  pivot_longer(cols = c(contains(c("_a"))),
               names_to = c("lithVisInfo", "lithSub"),
               names_sep = "_",
               values_to = 'lithValue') %>%
  filter(predVisInfo == lithVisInfo) %>%
  select(-c(lithVisInfo, lithSub)) %>%
  mutate(dataset = 'Average') %>% 
  rbind(.,meanCoordsChrom)

###* varying nearest % 
library(ggpmisc)
library(patchwork)

predVisInfoLabeller <- c('lumCoord' = 'Luminance',
                         xCoord = 'G:R opponent', yCoord = 'Y:B opponent')
statPoly <- stat_poly_eq(
  formula = y~x, 
  aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label..,
                    sep = "~~~")), parse = TRUE, vstep = 3, size =3)
substrateLabeller = c('a' = 'Lithops', 'b' = 'Rock', 'c' = 'Soil')

p1<-meanAll %>% 
  filter(predVisInfo == 'lumCoord') %>%
  ggplot(aes(x = predValue, y = lithValue))+
  geom_point()+
  geom_abline(lty='dashed')+
  geom_smooth(method='lm', colour = '#75263b', se = F, lwd = .7)+
  #coord_equal()+
  #xlab(NULL)+
  xlab('Background subs. mean luminance')+
  ylab('Lithops geometric mean luminance')+
  statPoly+ 
  facet_grid(dataset~predSub,
             labeller = labeller(predSub = substrateLabeller))+
  theme(strip.background = element_rect(fill = 'grey91'))

p2 <- meanAll %>% 
  filter(predVisInfo == 'xCoord') %>%
  ggplot(aes(x = predValue, y = lithValue))+
  geom_point()+
  geom_abline(lty='dashed')+
  geom_smooth(method='lm', colour = '#75263b', se = F, lwd = .7)+
  #coord_equal()+
  #xlab(NULL)+
  xlab('Background subs. mean G:R opponent')+
  ylab('Lithops mean G:R opponent')+
  statPoly+ 
  facet_grid(dataset~predSub,
             labeller = labeller(predSub = substrateLabeller))+
  theme(strip.background = element_rect(fill = 'grey91'))

p3 <- meanAll %>% 
  filter(predVisInfo == 'yCoord') %>%
  ggplot(aes(x = predValue, y = lithValue))+
  geom_point()+
  geom_abline(lty='dashed')+
  geom_smooth(method='lm', colour = '#75263b', se = F, lwd = .7)+
  #coord_equal()+
  #xlab(NULL)+
  xlab('Background subs. mean Y:B opponent')+
  ylab('Lithops mean Y:B opponent')+
  statPoly+ 
  facet_grid(dataset~predSub,
             labeller = labeller(predSub = substrateLabeller))+
  theme(strip.background = element_rect(fill = 'grey91'))

png("output\\rock-v-soil\\L-coords-af-RS-coords-amalg.png",
    type = 'cairo', units = "mm", res = 300,
    width = 215, height = 260)
p1+p2+p3 + plot_layout(nrow = 3)
dev.off()
rm(custFacet, custTheme,statPoly, gMeanLum,
   meanCoordsChrom, nrstLithops, p1,p2,p3)

# account for non indep of populations and species w random factors
library(lme4)
library(car)
meanAll <- merge(meanAll, unique(df[,c('abbrevs', 'collapseSpp')]))


###** Facet by spp 

meanAvg <- meanAll %>%
  mutate(facetVar = paste(predVisInfo,predSub)) %>%
  filter(dataset == "Average") %>%
  # cbind(.,props[match(.$abbrevs, props$abbrevs),
  #               c("rock_cover", "soil_cover")]) %>%
  cbind(., df[match(.$abbrevs, df$abbrevs),
              c("collapseSpp", "pop_spp")]) %>%
  separate(., abbrevs, "siteLabs", extra = 'drop')

# meanCoords <- x %>% pivot_longer(cols = contains(c("lithops")),
#                                  names_to = c("lithVisInfo", "lithSub"),
#                                  names_sep = "_",
#                                  values_to = 'lithValue') %>%
#   mutate(facetVar = paste(predVisInfo,predSub)) %>%
#   filter(predVisInfo == lithVisInfo, nrstSet == "100") %>%
#   select(-c(lithSub, lithVisInfo)) %>%
#   cbind(.,props[match(.$abbrevs, props$abbrevs),
#                 c("rock_cover", "soil_cover")]) %>%
#   cbind(., df[match(.$abbrevs, df$abbrevs),
#               c("collapseSpp", "pop_spp")]) %>%
#   separate(., abbrevs, "siteLabs", extra = 'drop')

library(ggforce)
library(cowplot)
library(ggrepel)

sppLabeller<- c(olivacea = 'L. olivacea', localis = 'L. localis',
                marmorataSpp = 'L. marmorata', divergensSpp = 'L. divergens',
                halliiSpp = 'L. hallii', leslieiSpp = 'L. lesliei',
                otzeniana = 'L. otzeniana', compt. = 'L. comptonii',
                fulleri = 'L. fulleri', bromf. = 'L. bromfieldii',
                hookeriSpp = 'L. hookeri', meyeri = 'L. meyeri',
                dint. = 'L. dinteri', verruc. = 'L. verruculosa',
                singlePops = 'singlePops')
substrateLabeller = c('a' = 'Lithops', 'b' = 'Rock', 'c' = 'Soil')

custText<- geom_text_repel(aes(label = siteLabs),
                           force = .5,
                           size = 3,
                           max.overlaps = Inf,
                           min.segment.length = 0.5)
custTheme <- theme(strip.placement = 'outside',
                   strip.background = element_blank(),
                   strip.text.y = element_blank(),
                   plot.title = element_text(hjust = 0.5),
                   legend.position = 'none')

p1 <- vector(mode = 'list', length = 2)
p2 <- vector(mode = 'list', length = 2)
p3 <- vector(mode = 'list', length = 2)

for(i in 1:2){
  j <- seq(1, length(unique(meanAvg$collapseSpp)), 10)

  pDataBg <- meanAvg %>%
    filter(dataset == 'Average',
           predVisInfo == 'lumCoord')
  pDataOver <- meanAvg[
    meanAvg$collapseSpp %in%
      unique(meanAvg$collapseSpp)[j[i]:(j[i]+9)],] %>%
    filter(dataset == 'Average',
           predVisInfo == 'lumCoord')

  p1[[i]]<- ggplot(pDataOver,aes(x= predValue, y=lithValue))+
    geom_point(data = pDataBg[, c("predSub","predValue", "lithValue")],
               color = "gray")+
    geom_point(size = 2, pch =21, fill = '#75263b')+
    # geom_point(aes(fill = rock_cover), size = 2, pch =21)+
    geom_abline(intercept = 0, lty = 2)+
    #labs(title = 'Luminance')+
    xlab('Substrate luminance')+
    ylab('Lithops luminance')+
    facet_grid(collapseSpp~predSub, switch = 'y',
               labeller = labeller(collapseSpp = sppLabeller,
                                   predSub = substrateLabeller))+
    custText+
    theme(strip.placement = 'outside',
          strip.text.y = element_text(face = 'italic'),
          strip.background = element_blank(),
          #panel.spacing.y = unit(0,"lines"),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none')

  pDataBg<-meanAvg %>%
    filter(dataset == 'Average',
           predVisInfo == 'xCoord')
  pDataOver <- meanAvg[
    meanAvg$collapseSpp %in%
      unique(meanAvg$collapseSpp)[j[i]:(j[i]+9)],] %>%
    filter(dataset == 'Average',
           predVisInfo == 'xCoord')

  p2[[i]]<- ggplot(pDataOver,aes(x= predValue, y=lithValue))+
    geom_point(data = pDataBg[, c('predSub',"predValue", "lithValue")],
               color = "gray")+
    geom_point(size = 2, pch =21, fill = '#75263b')+
    # geom_point(aes(fill = rock_cover), size = 2, pch =21)+
    geom_abline(intercept = 0, lty = 2)+
    custText+
    #labs(title = 'G:R opponent')+
    xlab('Substrate G:R opponent')+
    ylab('Lithops G:R opponent')+
    facet_grid(collapseSpp~predSub, 
               labeller = labeller(predSub = substrateLabeller))+
    custTheme

  pDataBg<- meanAvg %>%
    filter(dataset == 'Average',
           predVisInfo == 'yCoord')
  pDataOver <- meanAvg[
    meanAvg$collapseSpp %in%
      unique(meanAvg$collapseSpp)[j[i]:(j[i]+9)],] %>%
    filter(dataset == 'Average',
           predVisInfo == 'yCoord')

  p3[[i]]<- ggplot(pDataOver,aes(x= predValue, y=lithValue))+
    geom_point(data = pDataBg[, c('predSub',"predValue", "lithValue")],
               color = "gray")+
    geom_point(size = 2, pch =21, fill = '#75263b')+
    # geom_point(aes(fill = rock_cover), size = 2, pch =21)+
    geom_abline(intercept = 0, lty = 2)+
    custText+
    #labs(title = 'Y:B opponent')+
    xlab('Substrate Y:B opponent')+
    ylab('Lithops Y:B opponent')+
    facet_grid(collapseSpp~predSub,
               labeller = labeller(predSub = substrateLabeller))+
    custTheme

}


png("output\\rock-v-soil\\L-coords-af-RS-coords-by-spp-amalg-1.png",
    type = 'cairo', units = "mm", res = 300,
    width = 215, height = 300)
plot_grid(p1[[1]],p2[[1]],p3[[1]], ncol = 3)
dev.off()
png("output\\rock-v-soil\\L-coords-af-RS-coords-by-spp-amalg-2.png",
    type = 'cairo', units = "mm", res = 300,
    width = 215, height = 150)
plot_grid(p1[[2]],p2[[2]],p3[[2]], ncol = 3)
dev.off()

rm(p1,p2,p3, meanAll, meanAvg)


# Multiple regression - predicting Lith-sub dists ####
# Using surface area, RS distance, LS dist, variance 

# get LR, LS, RS distances

nrstSub_bc <- nrstSub_bc %>% filter(mspec.x == mspec.y)

nrstLithops <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                   "abbrevs", "mspec","roi", "substrate")]),
                     unique(nrstSub_bc[, c("abbrevs.x", "mspec.x", "roi.x")]),
                     by.y = c("abbrevs.x","mspec.x","roi.x"),
                     by.x = c("abbrevs","mspec", "roi")) %>%
  distinct(.) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean)

gMeanLum <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                "abbrevs", "mspec","roi", "substrate")]),
                  unique(nrstSub_bc[
                    nrstSub_bc$visInfo == 'dL',
                    c("abbrevs.x", "mspec.x", "roi.y")]),
                  by.y = c("abbrevs.x","mspec.x","roi.y"),
                  by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(abbrevs, substrate, mspec) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  rbind(nrstLithops,.) %>%
  as.data.frame(.) %>%
  split(.,.$abbrevs)

gMeanChrom <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                  "abbrevs", "mspec","roi", "substrate")]),
                    unique(nrstSub_bc[
                      nrstSub_bc$visInfo == 'dS',
                      c("abbrevs.x", "mspec.x", "roi.y")]),
                    by.y = c("abbrevs.x","mspec.x","roi.y"),
                    by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(abbrevs, substrate, mspec) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  rbind(nrstLithops,.) %>%
  as.data.frame(.) %>%
  split(.,.$abbrevs)

gMeanAll <- df %>%
  group_by(abbrevs, substrate, mspec) %>% # each image contributes equally to gmean dist
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  as.data.frame(.) %>%
  split(.,.$abbrevs)

gMeanLum <- lapply(gMeanLum, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

gMeanChrom <- lapply(gMeanChrom, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

gMeanAll <- lapply(gMeanAll, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

gMeanLum <- lapply(gMeanLum, coldist,
                   n = c(1,16,32), weber = .05, weber.achro = .05,
                   weber.ref = "longest",
                   qcatch = "Qi", achromatic = TRUE) %>%
  rbindlist(., idcol = 'abbrevs') %>%
  rename(chroma = dS, lum = dL) %>%
  pivot_longer(cols = c(chroma,lum),
               names_to = 'visInfo', values_to = 'distance') %>%
  filter(patch1 == 'a' | patch2 == 'a', visInfo == 'lum') %>%
  mutate(patch1 = fct_recode(patch1, L = 'a', R = 'b', S = 'c'),
         patch2 = fct_recode(patch2, L = 'a', R = 'b', S = 'c'),
         contrast = paste(patch1, patch2,sep = '')) %>%
  select(-c(patch2,patch1)) %>%
  mutate(dataset = 'nrst')

gMeanChrom <- lapply(gMeanChrom, coldist,
                     n = c(1,16,32), weber = .05, weber.achro = .05,
                     weber.ref = "longest",
                     qcatch = "Qi", achromatic = TRUE) %>%
  rbindlist(., idcol = 'abbrevs') %>%
  rename(chroma = dS, 'lum' = dL) %>%
  pivot_longer(cols = c(chroma,lum),
               names_to = 'visInfo', values_to = 'distance') %>%
  filter(patch1 == 'a' | patch2 == 'a', visInfo == 'chroma') %>%
  mutate(patch1 = fct_recode(patch1, L = 'a', R = 'b', S = 'c'),
         patch2 = fct_recode(patch2, L = 'a', R = 'b', S = 'c'),
         contrast = paste(patch1, patch2,sep = '')) %>%
  select(-c(patch2,patch1)) %>%
  mutate(dataset = 'nrst')

gMeanAll <- lapply(gMeanAll, coldist,
                   n = c(1,16,32), weber = .05, weber.achro = .05, 
                   weber.ref = "longest",
                   qcatch = "Qi", achromatic = TRUE) %>%
  rbindlist(., idcol = 'abbrevs') %>%
  mutate(patch1 = fct_recode(patch1, L = 'a', R = 'b', S = 'c'),
         patch2 = fct_recode(patch2, L = 'a', R = 'b', S = 'c'),
         contrast = paste(patch1, patch2,sep = '')) %>%
  select(-c(patch1,patch2)) %>%
  mutate(dataset = 'avg') %>%
  rename(.,chroma = dS, 'lum' = dL) %>%
  pivot_longer(cols = c('chroma', 'lum'),
               values_to = 'distance',
               names_to = 'visInfo')

gMeans <- rbind(gMeanAll, gMeanLum, gMeanChrom) %>%
  pivot_wider(names_from = 'contrast', values_from = 'distance')
rm(gMeanAll, gMeanLum, gMeanChrom, nrstLithops)

# Get substrate freq.

# Add substrate frequency variable
source("rock-v-soil//rock-estimate.R")
props <- dfRockEst %>% count(abbrevs, substrate) %>%
  pivot_wider(names_from = substrate, values_from = n) %>%
  mutate(rockCover = b/(b+c)*100,
         soilCover = c/(b+c)*100) %>%
  select(-c(a,b,c))
source("loading-cleaning.R") # not ideal. could change name of rock estimate df

# get sd
visHet <- df %>%
  group_by(abbrevs,substrate) %>%
  mutate(lumCoord = log(lumMean)/0.05) %>%
  summarise_at(vars(lumCoord, xCoord, yCoord), sd) %>%
  rename(lumCoordSd = 'lumCoord', xCoordSd = 'xCoord',
         yCoordSd = 'yCoord') %>%
  mutate(substrate = fct_recode(substrate, L = 'a', R = 'b', S = 'c'),
         chromaSd = (xCoordSd+yCoordSd)/2) %>%
  pivot_wider(names_from = substrate,
              values_from = c(lumCoordSd, xCoordSd, yCoordSd,chromaSd)) %>%
  select(-c(contains('_L')))

# merge datasets
glmDf <- merge(visHet,props) %>%
  merge(., gMeans) %>%
  pivot_wider(names_from = c(dataset,visInfo),
              values_from = c(LR,LS,RS),
              names_glue = "{dataset}_{visInfo}_{.value}") %>%
  mutate(avg_chroma_LRLS = avg_chroma_LR - avg_chroma_LS,
         avg_lum_LRLS = abs(avg_lum_LR - avg_chroma_LS))

# Multiple regressions

avg_lum <- lm(data = glmDf,
              avg_lum_LR~avg_lum_LS+avg_lum_RS+rockCover+lumCoordSd_R)
avg_chroma <- lm(data = glmDf,
                 avg_chroma_LR~avg_chroma_LS+
                   avg_chroma_RS+rockCover+chromaSd_R)
nrst_lum <- lm(data = glmDf,
               nrst_lum_LR~ nrst_lum_LS+ avg_lum_RS+
                 rockCover+ lumCoordSd_R)
nrst_chroma <- lm(data = glmDf, 
                  nrst_chroma_LR~ nrst_chroma_LS+ avg_chroma_RS+
                    rockCover+ chromaSd_R)
rel_index_lum <- lm(data = glmDf,
                    avg_lum_LRLS~avg_lum_LS+avg_lum_RS+
                      rockCover+lumCoordSd_R)
rel_index_chroma <- lm(data = glmDf,
                    avg_chroma_LRLS~avg_chroma_LS+avg_chroma_RS+
                      rockCover+chromaSd_R)

anova(avg_chroma)
anova(nrst_chroma)
anova(avg_lum)
anova(nrst_lum)
summary(avg_chroma)
summary(avg_lum)
summary(nrst_chroma)
summary(nrst_lum)
summary(rel_index_lum)
summary(rel_index_chroma)

library(stargazer)
stargazer(avg_lum,nrst_lum,avg_chroma,nrst_chroma,
          out = 'output//rock-v-soil//tables//predicting-RS-Reg-table.htm',
          star.char = c("*", "**", "***"),
          star.cutoffs = c(.05, .01, .001))

rm(nrst_chroma, nrst_lum,avg_chroma,avg_lum,
   rel_index_chroma, rel_index_lum,visHet,props,gMeans,
   dfRockEst)
# Plots
# cover predicts distance
# statPoly <- stat_poly_eq(
#   formula = y~x,aes(label = paste(..eq.label..,..rr.label..,
#                                   ..p.value.label.., sep = "~~~~")),
#                            parse = TRUE,
#   size = 2)
# geomSmoothCust <- geom_smooth(method = "lm", color = "black",
#                                 formula = y~x, se =F, lwd = .5)
library(egg)

LRlumPredict <- glmDf %>%
  select(c(avg_lum_LR,avg_lum_LS, avg_lum_RS, rockCover, lumCoordSd_R)) %>%
  pivot_longer(cols = c(avg_lum_LS, avg_lum_RS, rockCover, lumCoordSd_R),
               names_to = 'LRlumDistPreds')

sigPredictors <- LRlumPredict %>%
  filter(LRlumDistPreds %in% c('avg_lum_LS',  'avg_lum_RS'))
  
p1<- ggplot(LRlumPredict, aes(value, avg_lum_LR))+
  geom_point()+
  geom_smooth(data = sigPredictors, 
              method = "lm", color = "#75263b",
              formula = y~x, se =F, lwd = .5)+
  facet_wrap(vars(LRlumDistPreds),
             scales = 'free_x', nrow = 2,
             strip.position = 'bottom',
             labeller = labeller(
               LRlumDistPreds = c(avg_lum_LS = 'LS luminance distance',
                                  avg_lum_RS = 'RS luminance distance',
                                  rockCover = 'Rock cover',
                                  lumCoordSd_R = 'Rock luminance sd')))+
   labs(x = '',
        y = 'LR luminance distance', 
        title = 'Predictors of Lithops-rock luminance distance')+
  theme(axis.title.y = element_text(size = 9),
        strip.background.x = element_blank(),
        plot.title = element_text(hjust = .5))
#p1 <- tag_facet(p1)+theme(strip.text = element_text())

rm(LRlumPredict, sigPredictors)
 
LRPredChroma <- glmDf %>%
   select(c(avg_chroma_LS,avg_chroma_LR, avg_chroma_RS,
            rockCover, chromaSd_R)) %>%
   pivot_longer(cols = c(avg_chroma_LS, avg_chroma_RS,
                         rockCover, chromaSd_R),
                names_to = 'chromaPreds') %>%
  mutate(chromaPreds =
           fct_relevel(chromaPreds, 'avg_chroma_LS', 'avg_chroma_RS',
                                  'chromaSd_R', 'rockCover'))

chromaSigPreds <- LRPredChroma %>%
  filter(chromaPreds %in% c('avg_chroma_RS','avg_chroma_LS', 'chromaSd_R'))

p2<- ggplot(LRPredChroma, aes(value,avg_chroma_LR))+
   geom_point()+
   facet_wrap(vars(chromaPreds),
              scales = 'free_x', nrow = 2,
              strip.position = 'bottom',
              labeller = labeller(
                chromaPreds = c(avg_chroma_LS = 'LS chroma distance',
                                   avg_chroma_RS = 'RS chroma distance',
                                   chromaVar_R = 'Rock chroma sd',
                                   rockCover = 'Rock cover')))+
   labs(x = '',
        y = 'LR chroma distance', 
        title = 'Predictors of Lithops-rock chroma distance')+ 
   theme(axis.title.y = element_text(size = 9),
         strip.background.x = element_blank(),
         plot.title = element_text(hjust = .5))+
   geom_smooth(data = chromaSigPreds, 
               method = "lm", color = "#75263b",
               formula = y~x, se =F, lwd = .5)

png('output//rock-v-soil//substrate-distance-predictors.png',
    width = 215,height = 215, res = 300, type = 'cairo', units = 'mm')
p1/p2+plot_layout()+plot_annotation(tag_levels = 'A')
dev.off()

# Specialising - distance~substrate*species + substrate*geology ####

nrstSub_bc <- nrstSub_bc %>% filter(mspec.x == mspec.y)

nrstLithops <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                   "abbrevs", "mspec","roi", "substrate")]),
                     unique(nrstSub_bc[, c("abbrevs.x", "mspec.x", "roi.x")]),
                     by.y = c("abbrevs.x","mspec.x","roi.x"),
                     by.x = c("abbrevs","mspec", "roi")) %>%
  distinct(.) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean)

gMeanLum <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                "abbrevs", "mspec","roi", "substrate")]),
                  unique(nrstSub_bc[
                    nrstSub_bc$visInfo == 'dL',
                    c("abbrevs.x", "mspec.x", "roi.y")]),
                  by.y = c("abbrevs.x","mspec.x","roi.y"),
                  by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(abbrevs, substrate, mspec) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  rbind(nrstLithops,.) %>%
  as.data.frame(.) %>%
  split(.,.$abbrevs)

gMeanChrom <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                  "abbrevs", "mspec","roi", "substrate")]),
                    unique(nrstSub_bc[
                      nrstSub_bc$visInfo == 'dS',
                      c("abbrevs.x", "mspec.x", "roi.y")]),
                    by.y = c("abbrevs.x","mspec.x","roi.y"),
                    by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(abbrevs, substrate, mspec) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  rbind(nrstLithops,.) %>%
  as.data.frame(.) %>%
  split(.,.$abbrevs)

gMeanAll <- df %>%
  group_by(abbrevs, substrate, mspec) %>% # each image contributes equally to gmean dist
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  as.data.frame(.) %>%
  split(.,.$abbrevs)

gMeanLum <- lapply(gMeanLum, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

gMeanChrom <- lapply(gMeanChrom, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

gMeanAll <- lapply(gMeanAll, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

gMeanLum <- lapply(gMeanLum, coldist,
                   n = c(1,16,32), weber = .05, weber.achro = .05,
                   weber.ref = "longest",
                   qcatch = "Qi", achromatic = TRUE) %>%
  rbindlist(., idcol = 'abbrevs') %>%
  rename(chroma = dS, lum = dL) %>%
  pivot_longer(cols = c(chroma,lum),
               names_to = 'visInfo', values_to = 'distance') %>%
  filter(patch1 == 'a' | patch2 == 'a', visInfo == 'lum') %>%
  mutate(patch1 = fct_recode(patch1, L = 'a', R = 'b', S = 'c'),
         patch2 = fct_recode(patch2, L = 'a', R = 'b', S = 'c'),
         contrast = paste(patch1, patch2,sep = '')) %>%
  select(-c(patch2,patch1)) %>%
  mutate(dataset = 'nrst')

gMeanChrom <- lapply(gMeanChrom, coldist,
                     n = c(1,16,32), weber = .05, weber.achro = .05,
                     weber.ref = "longest",
                     qcatch = "Qi", achromatic = TRUE) %>%
  rbindlist(., idcol = 'abbrevs') %>%
  rename(chroma = dS, 'lum' = dL) %>%
  pivot_longer(cols = c(chroma,lum),
               names_to = 'visInfo', values_to = 'distance') %>%
  filter(patch1 == 'a' | patch2 == 'a', visInfo == 'chroma') %>%
  mutate(patch1 = fct_recode(patch1, L = 'a', R = 'b', S = 'c'),
         patch2 = fct_recode(patch2, L = 'a', R = 'b', S = 'c'),
         contrast = paste(patch1, patch2,sep = '')) %>%
  select(-c(patch2,patch1)) %>%
  mutate(dataset = 'nrst')

gMeanAll <- lapply(gMeanAll, coldist,
                   n = c(1,16,32), weber = .05, weber.achro = .05, 
                   weber.ref = "longest",
                   qcatch = "Qi", achromatic = TRUE) %>%
  rbindlist(., idcol = 'abbrevs') %>%
  mutate(patch1 = fct_recode(patch1, L = 'a', R = 'b', S = 'c'),
         patch2 = fct_recode(patch2, L = 'a', R = 'b', S = 'c'),
         contrast = paste(patch1, patch2,sep = '')) %>%
  select(-c(patch1,patch2)) %>%
  mutate(dataset = 'avg') %>%
  rename(.,chroma = dS, 'lum' = dL) %>%
  pivot_longer(cols = c('chroma', 'lum'),
               values_to = 'distance',
               names_to = 'visInfo') 

gMeans <- rbind(gMeanAll, gMeanLum, gMeanChrom) %>%
  filter(contrast != 'RS') %>%
  merge(.,unique(df[,c('abbrevs', 'collapseSpp', 'allan_geol')])) %>%
  pivot_wider(values_from = 'distance', names_from = c('visInfo'))

lumAvgMod <- lm(data = subset(gMeans, dataset == 'avg'),
                      lum~contrast*collapseSpp)
lumNrstMod <- lm(data = subset(gMeans, dataset == 'nrst'),
                 lum~contrast*collapseSpp)
chromaAvgMod <- lm(data = subset(gMeans, dataset == 'avg'),
                   chroma~contrast*collapseSpp)
chromaNrstMod <- lm(data = subset(gMeans, dataset == 'nrst'),
                    chroma~contrast*collapseSpp)

# Coefficients
stargazer(lumAvgMod, lumNrstMod, chromaAvgMod, chromaNrstMod,
          single.row = T,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digits = 2,
          out = 'output//rock-v-soil//tables//table-interac-sub&spp-dist-betas.htm')

# Anova tables
anova_results <- lapply(list(lumAvgMod, lumNrstMod,
                             chromaAvgMod, chromaNrstMod), anova) %>%
  rbindlist() %>%
  cbind( c(rep(c('Substrate', 'Species',
                 'Substrate:species', 'Residuals'), 4)),.)
colnames(anova_results) <- c("Sum Sq", "Df", "F value", "Pr(>F)","") 
row.names(anova_results) <- NULL

library(kableExtra)
anova_results %>% kable("html", digits=3) %>% 
  kable_styling(bootstrap_options = "basic", full_width = F) %>% 
  pack_rows(., "Luminance distance (avg)", 1, 4) %>%
  pack_rows(., "Luminance distance (nrst)", 5, 8) %>%
  pack_rows(., "Chroma distance (avg)", 9, 12) %>%
  pack_rows(., "Chroma distance (nrst)", 13, 16) %>%
  save_kable('output//rock-v-soil//tables//table-interac-spp&sub-pred-dist.htm')

rm(gMeans, anova_results, lumAvgMod, lumNrstMod, 
   chromaAvgMod, chromaNrstMod)

# Phylogenetic signal - ANOVAs on Lithops coordinates ####

lithopsCoords <- df %>%
  filter(substrate == 'a') %>%
  group_by(abbrevs) %>%
  summarise(across(contains("Coord"), mean, .names = "mean_{.col}")) %>%
  merge(., unique(df[,c('abbrevs', 'collapseSpp')]))

phyloSignalLum <- lm(data = lithopsCoords, mean_lumCoord~collapseSpp)
phyloSignalxCoord <- lm(data = lithopsCoords, mean_xCoord~collapseSpp)
phyloSignalyCoord <- lm(data = lithopsCoords, mean_yCoord~collapseSpp)

aovTabxCoord<- anova(phyloSignalxCoord)
aovTabyCoord<- anova(phyloSignalyCoord)
aovTabLumCoord<- anova(phyloSignalLum)

summary(phyloSignalLum)
summary(phyloSignalxCoord)
summary(phyloSignalyCoord)

anova_results <- data.frame(
  cbind(c("Species", "Residuals",
          "Species", "Residuals",
          "Species", "Residuals"), 
        rbind(aovTabLumCoord, aovTabxCoord, aovTabyCoord))) 
colnames(anova_results) <- c("", "Df",  "Sum Sq", "Mean Sq",
                             "F value", "Pr(>F)")
row.names(anova_results) <- NULL

library(kableExtra)
anova_results %>% kable("html", digits=3) %>% 
  kable_styling(bootstrap_options = "basic", full_width = F) %>% 
  pack_rows(., "Lithops luminance", 1, 2) %>%
  pack_rows(., "Lithops G:R opponent", 3, 4) %>%
  pack_rows(., "Lithops Y:B opponent", 5, 6) %>%
  save_kable('output//rock-v-soil//tables//anova-table-LithopsCoord-phyloSignal.htm')

rm(lithopsCoords, phyloSignalLum)


# Background visual properties ####

# ttest rock sd versus soil sd

bg_propties <- df %>%
  group_by(abbrevs, substrate) %>%
  summarise(across(c(xCoord, yCoord, lumCoord),
                   c(mean = mean, sd = sd),
                   .names = "{.fn}_{.col}")) 

t_test <- function(df, alt = 'two.sided', paired = T){
  tidy(t.test(df$b, df$c))
}

rvsVar <- bg_propties %>%
  filter(substrate != 'a') %>%
  select(c(abbrevs,substrate,contains('sd'))) %>%
  pivot_longer(cols = c(where(is.numeric)),
               names_to = c('statistic', 'visInfo'),
               names_sep = '_',
               values_to = 'sd') 

pairedTests <- rvsVar %>%
  pivot_wider(names_from = 'substrate', 
              values_from = 'sd') %>%
  group_by(visInfo) %>%
  nest() %>%
  mutate(ttest = map(data, t_test)) %>%
  unnest(ttest) %>%
  mutate(labelz = 
           case_when(p.value < 0.001 ~ '***',
                     p.value < 0.01 ~ '**',
                     p.value < 0.05 ~ '*'))

p7<- ggplot(rvsVar, 
         aes(y = sd, x = substrate))+
  geom_boxplot(aes(fill = substrate))+
  geom_text(data = pairedTests,
            aes(label = labelz, x = 1.5, y = Inf), vjust = 1.2, size = 8)+
  facet_wrap(vars(visInfo), scales = 'free', ncol = 3,
             labeller = labeller(visInfo = c(
               xCoord = 'G:R opponent',
               yCoord = 'Y:B opponent',
               lumCoord = 'Luminance')))+
  scale_fill_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"))+
  scale_x_discrete(labels = c('Rock', 'Soil'))+
  theme(legend.position = 'none',
        strip.background = element_rect(fill = 'grey91'))
png('output//rock-v-soil//RS-var-relship-boxplots.png',
     width = 215, height = 100, res = 300, type = 'cairo', units = 'mm')
p7
dev.off()

rm(rvsVar, pairedTests)

# Plot
bg_propties <- bg_propties %>%
  pivot_wider(names_from = substrate,
              values_from = c(where(is.numeric)))

statPoly <- stat_poly_eq(
  formula = y~x,
  aes(label = paste(..eq.label..,
                    ..rr.label..,
                    ..p.value.label..,
                    sep = "~~~~")), parse = TRUE, size = 2)

p1<-ggplot(bg_propties, aes(sd_lumCoord_b, sd_lumCoord_c))+
  geom_point()+
  coord_equal(ratio = 1)+
  geom_abline(lty=2)+
  #  stat_summary(fun.data= mean_cl_normal)+ 
  #geom_smooth(method='lm', colour = '#75263b', se = F, lwd = .5)+
  scale_x_continuous(
    limits = c(min(bg_propties$sd_lumCoord_b,
                   bg_propties$sd_lumCoord_c),
                    max(bg_propties$sd_lumCoord_b,
                        bg_propties$sd_lumCoord_c)))+
  scale_y_continuous(
    limits = c(min(bg_propties$sd_lumCoord_b,
                   bg_propties$sd_lumCoord_c),
               max(bg_propties$sd_lumCoord_b,
                   bg_propties$sd_lumCoord_c)))+
  scale_color_viridis_c(option = 'magma')+
  xlab("Rock luminance sd")+
  ylab("Soil luminance sd")+
  statPoly+
  theme(legend.position = 'none')

p2<-ggplot(bg_propties, aes(sd_xCoord_b, sd_xCoord_c))+
  geom_point()+
  coord_equal(ratio = 1)+
  geom_abline(lty=2)+
  #  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = '#75263b', se = F, lwd = .5)+
  scale_x_continuous(
      limits = c(min(bg_propties$sd_xCoord_b,
                     bg_propties$sd_xCoord_c),
                 max(bg_propties$sd_xCoord_b,
                     bg_propties$sd_xCoord_c)))+
      scale_y_continuous(
        limits = c(min(bg_propties$sd_xCoord_b,
                       bg_propties$sd_xCoord_c),
                   max(bg_propties$sd_xCoord_b,
                       bg_propties$sd_xCoord_c)))+
  scale_color_viridis_c(option = 'magma')+
  xlab("Rock G:R sd")+
  ylab("Soil G:R sd")+
  statPoly+
  theme(legend.position = 'none')

p3<-ggplot(bg_propties, aes(sd_yCoord_b, sd_yCoord_c))+
  geom_point()+
  coord_equal(ratio = 1)+
  geom_abline(lty=2)+
  #  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = '#75263b', se = F, lwd = .5)+
  scale_x_continuous(
    limits = c(min(bg_propties$sd_yCoord_b,
                   bg_propties$sd_yCoord_c),
               max(bg_propties$sd_yCoord_b,
                   bg_propties$sd_yCoord_c)))+
  scale_y_continuous(
    limits = c(min(bg_propties$sd_yCoord_b,
                   bg_propties$sd_yCoord_c),
               max(bg_propties$sd_yCoord_b,
                   bg_propties$sd_yCoord_c)))+
  scale_color_viridis_c(option = 'magma')+
  xlab("Rock Y:B sd")+
  ylab("Soil Y:B sd")+
  statPoly+
  theme(legend.position = 'none')

# Mean coordinate relship plots
rvsMeanCoords <-  df %>% 
  filter(substrate != 'a') %>%
  mutate(lumMeanCoord = log(lumMean)/0.05) %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(xCoord,yCoord,lumMeanCoord), mean) %>%
  pivot_wider(names_from = substrate,
              values_from = c('xCoord', 'yCoord', 'lumMeanCoord'))
  # pivot_longer(cols = c(xCoord,yCoord,lumMeanCoord),
  #              names_to = 'visInfo', values_to = 'coordinate')

p4<-ggplot(bg_propties, aes(mean_lumCoord_b, mean_lumCoord_c))+
  geom_point()+
  coord_equal(ratio = 1)+
  geom_abline(lty=2)+
  #stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = '#75263b', se = F, lwd = .5)+
  scale_x_continuous(
    limits = c(min(bg_propties$mean_lumCoord_b,
                   bg_propties$mean_lumCoord_c),
               max(bg_propties$mean_lumCoord_b,
                   bg_propties$mean_lumCoord_c)))+
  scale_y_continuous(
    limits = c(min(bg_propties$mean_lumCoord_b,
                   bg_propties$mean_lumCoord_c),
               max(bg_propties$mean_lumCoord_b,
                   bg_propties$mean_lumCoord_c)))+
  scale_color_viridis_c(option = 'magma')+
  xlab("Rock luminance")+
  ylab("Soil luminance")+
  statPoly
  #theme(legend.position = 'none')

p5<-ggplot(bg_propties, aes(mean_xCoord_b, mean_xCoord_c))+
  geom_point()+
  coord_equal(ratio = 1)+
  geom_abline(lty=2)+
  #stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = '#75263b', se = F, lwd = .5)+
  scale_x_continuous(
    limits = c(min(bg_propties$mean_xCoord_b,
                   bg_propties$mean_xCoord_c),
               max(bg_propties$mean_xCoord_b,
                   bg_propties$mean_xCoord_c)))+
  scale_y_continuous(
    limits = c(min(bg_propties$mean_xCoord_b,
                   bg_propties$mean_xCoord_c),
               max(bg_propties$mean_xCoord_b,
                   bg_propties$mean_xCoord_c)))+
  scale_color_viridis_c(option = 'magma')+
  xlab("Rock G:R opponent")+
  ylab("Soil G:R opponent")+
  statPoly+
  theme(legend.position = 'none')

p6<-ggplot(bg_propties, aes(mean_yCoord_b, mean_yCoord_c))+
  geom_point()+
  coord_equal(ratio = 1)+
  geom_abline(lty=2)+
  #stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = '#75263b', se = F, lwd = .5)+
  scale_x_continuous(
    limits = c(min(bg_propties$mean_yCoord_b,
                   bg_propties$mean_yCoord_c),
               max(bg_propties$mean_yCoord_b,
                   bg_propties$mean_yCoord_c)))+
  scale_y_continuous(
    limits = c(min(bg_propties$mean_yCoord_b,
                   bg_propties$mean_yCoord_c),
               max(bg_propties$mean_yCoord_b,
                   bg_propties$mean_yCoord_c)))+
  scale_color_viridis_c(option = 'magma')+
  xlab("Rock Y:B opponent")+
  ylab("Soil Y:B opponent")+
  statPoly+
  theme(legend.position = 'none')

png('output//rock-v-soil//rock-soil-coord&var-relship.png',
    width = 215, height = 160, res = 300, type = 'cairo', units = 'mm')
(p4+p5+p6+p1+p2+p3)+plot_layout()+ 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 8))
dev.off()

rm(bg_propties, p1,p2,p3,p4,p5,p6,p7,pairedTests)
# Within-pop Background differences (RvS) ####
# nearest rock to soil will be a bit more effort
logQi <- df %>% 
  group_by(abbrevs, mspec, substrate) %>%
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) %>%
  mutate(splitFct = paste(abbrevs,mspec, sep = '--')) %>%
  as.data.frame() %>% ungroup() %>%
  split(., .$splitFct)

logQi <- lapply(logQi, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

logQi <- Filter(function(x) nrow(x)>=2, logQi)
df %>%
  group_by(abbrevs, mspec) %>%
  filter(n_distinct(substrate)<2) 

perImageDistsAvg <- lapply(logQi, coldist,
                           n = c(1,16,32), weber = .05, weber.achro = .05, 
                           weber.ref = "longest",
                           qcatch = "Qi", achromatic = TRUE)

perImageDistsAvg <- rbindlist(perImageDistsAvg, idcol= "names") %>%
  separate(., names, into = c("abbrevs", "mspec"), sep = "--") %>%
  filter(patch1 != 'a' & patch2 != 'a') %>%
  merge(., unique(df[,c('abbrevs', 'collapseSpp')]),
        by.x = 'abbrevs', by.y = 'abbrevs') %>%
  pivot_longer(cols = c(dS,dL), 
               names_to = 'visInfo', values_to = 'distance') 

t_test <- function(df, mu = 1, alt = 'two.sided', conf.level = .95){
  tidy(t.test(df$distance,
         mu = mu,
         alt = alt,
         conf.level = conf.level))
}

ttestRVS <- perImageDistsAvg %>%
  group_by(abbrevs, visInfo) %>%
  nest() %>%
  mutate(ttest = map(data, t_test)) %>%
  unnest(ttest) %>%
  group_by(visInfo) %>%
  mutate(p.adj = p.adjust(p.value, n = 56, method = 'holm')) %>%
  merge(., unique(df[, c('abbrevs', 'collapseSpp')])) %>%
  mutate(labelPadj = case_when(p.adj < 0.05 ~ "*"),
         labelPval = case_when(p.value < 0.05 ~ '*'))
# Plot

# 
# facetVarLabeller <- c('dL Nearest 10%' = 'Lum Nearest 10%',
#                       'dS Nearest 10%' = 'Chroma Nearest 10%',
#                       'dL Average' = 'Luminance Average',
#                       'dS Average' = 'Chroma Average')
sppLabeller<- c(olivacea = 'L. oliv.', localis = 'L. loc.',
                marmorataSpp = 'L. marmorata', divergensSpp = 'L. divergens',
                halliiSpp = 'L. hallii', leslieiSpp = 'L. lesliei',
                otzeniana = 'L. otz.', compt. = 'L. compt.',
                fulleri = 'L. fulleri', bromf. = 'L. bromf.',
                hookeriSpp = 'L. hookeri', meyeri = 'L. me.', 
                dint. = 'L. din.', verruc. = 'L. verru.', 
                singlePops = 'singlePops')

visInfoRename = c(dL = 'RS avg lum dist', dS = 'RS avg chroma dist')

png(file = "output//rock-v-soil//RvS-sigDiff.png", type = 'cairo',
    width = 215, height = 270, units = "mm", res = 300) 
ggplot(ttestRVS,
       aes(y=estimate, x=abbrevs))+
  geom_point()+
  geom_errorbar(aes(x=abbrevs, ymax=conf.high, ymin=conf.low),
                position = position_dodge(width = .5))+
  geom_text(aes(y=Inf, label=labelPadj),
            colour = '#75263b',hjust = 1.5, vjust = 0.8)+
  geom_hline(yintercept = 1, colour = '#75263b', lty = 2)+
  coord_flip()+
  facet_grid(collapseSpp~visInfo, 
             switch = "y", scales = "free", space = "free_y",
             labeller = labeller(visInfo = visInfoRename,
                                 collapseSpp = sppLabeller))+
  ylab('Noise scaled distance (JNDs)')+
  theme(strip.placement = 'outside',
        strip.text.y = element_text(face = 'italic'),
        strip.background = element_rect(fill = 'grey91'),
        panel.spacing.y = unit(0,'lines'),
        legend.position = 'top')
dev.off()

rm(ttestRVS)
# bootstraps ####

bootPairedCoords <- df %>% as.data.frame(.) %>%  
  group_by(abbrevs, mspec, substrate) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), mean) %>% # ORDER NBNBNB!!!
  as.data.frame(.)
rownames(bootPairedCoords) <- make.unique(bootPairedCoords$substrate)
lsAbbrevs <- split(bootPairedCoords,bootPairedCoords$abbrevs)

lsSubVec <- lapply(lsAbbrevs, function(x){
  x <- x[, "substrate"]
}) 
lsAbbrevs <- lapply(lsAbbrevs, function(x) {
  x <- x[,which(names(x) %in% c("swMean", "mwMean", "lwMean", "lumMean"))] ; x
})
lsAbbrevs <- lapply(lsAbbrevs, function(x) {
  colnames(x)[which(names(x) %in% c("swMean", "mwMean", "lwMean", "lumMean"))] <- 
    c("s", "m", "l", "lum") ; x
})

bootRes <- vector(mode = "list", length = length(lsAbbrevs))
bootRes <- lapply(seq_along(lsAbbrevs),
                  function(x) bootcoldist(lsAbbrevs[[x]], 
                                          by = lsSubVec[[x]],
                                          n = c(1,16,32),
                                          weber = .05,
                                          weber.achro = .05, 
                                          weber.ref = "longest",
                                          qcatch = "Qi",
                                          achromatic = TRUE))
names(bootRes) <- names(lsAbbrevs)
bootResLs<- lapply(bootRes, function(x){
  x <- as.data.frame(x)
  x <- rownames_to_column(x, var = 'contrast')
})

bootResDf <- rbindlist(bootResLs, idcol = 'abbrevs') %>%
  mutate(contrast = fct_recode(.$contrast,
                               'LR' = "a-b",
                               'LS' = 'a-c',
                               'RS' = 'b-c')) %>% 
  pivot_longer(cols = contains(c('mean', 'lwr', 'upr')),
               names_to = c('visInfo','.value'),
               names_sep = '\\.') %>%
  cbind(., df[match(.$abbrevs, df$abbrevs),
              c("collapseSpp", "pop_spp")]) %>%
  arrange(desc(collapseSpp)) %>%
  mutate(abbrevs = as_factor(abbrevs))
levels()
# group_by(abbrevs) %>%
# mutate(incrCertainty = lwr[contrast == 'LR' & visInfo == 'dL']-
#          upr[contrast == 'LS' & visInfo == 'dS'],
#        abbrevs = fct_reorder(abbrevs, incrCertainty)) %>% 
# ungroup %>% mutate(abbrevs = fct_reorder(abbrevs, incrCertainty))

# Check if 95% conf intervals overlap
sigDifdS <- bootResDf %>%
  pivot_wider(names_from = c(contrast, visInfo),
              values_from = c('mean', 'lwr','upr'))  %>%
  filter(upr_LR_dS < lwr_LS_dS | upr_LS_dS < lwr_LR_dS) %>%
  pivot_longer(cols = contains(c('lwr','upr','mean')),
               names_to = c('.value', 'contrast', 'visInfo'),
               names_sep = '_') %>% 
  filter(visInfo != 'dL') %>%
  mutate(labelPos = 2.5)

annoDF <- bootResDf %>% 
  pivot_wider(names_from = c(contrast, visInfo),
              values_from = c('mean', 'lwr','upr'))  %>%
  filter(upr_LR_dL < lwr_LS_dL | upr_LS_dL < lwr_LR_dL) %>%
  pivot_longer(cols = contains(c('lwr','upr','mean')),
               names_to = c('.value', 'contrast', 'visInfo'),
               names_sep = '_') %>% 
  filter(visInfo != 'dS') %>%
  mutate(labelPos = 15) %>%
  rbind(., sigDifdS) %>%
  mutate(label = '*')

# plot
visInfoRename = c(dL = 'Lum dist', dS = 'Chroma dist')
png(file = "output//rock-v-soil//RvS-bootstraps-All.png", type = 'cairo',
    width = 200, height = 300, units = "mm", res = 300) 
ggplot(bootResDf %>% filter(contrast != 'RS' ),
       aes(y=mean, x=abbrevs,
           group=contrast,
           colour=contrast,
           fill=contrast))+
  geom_point(pch = 21, size = 1,
             position = position_dodge(width = .5))+
  geom_errorbar(aes(x=abbrevs, ymax=upr, ymin=lwr),
                position = position_dodge(width = .5))+
  geom_text(data=annoDF %>% filter(contrast != 'RS'),
            aes(y=labelPos, label=label), colour = 'black')+
  coord_flip()+
  facet_wrap(vars(visInfo), scales = 'free_x',
             strip.position = 'top',
             labeller = labeller(visInfo = visInfoRename))
dev.off()

# summary stats
annoDF %>% pivot_wider(names_from = c('contrast'),
                       values_from = c('mean', 'upr', 'lwr')) %>%
  filter(visInfo == 'dS' & mean_LR > mean_LS) 
annoDF %>% pivot_wider(names_from = c('contrast'),
                       values_from = c('mean', 'upr', 'lwr')) %>%
  filter(visInfo == 'dS' & mean_LR < mean_LS) 

bootResDf %>% pivot_wider(names_from = c(visInfo, contrast),
                          values_from = contains(c('lwr','upr','mean'))) %>%
  filter(mean_dL_LR > mean_dL_LS)

rm(list=ls(pattern="^boot"), annoDF, lsAbbrevs, lsSubVec, sigDifdS)

rm(distSpreadAll, p1,p2,p3)



# Other stuff etc. ####

# outliers
outliers <- jnds_dec[jnds_dec$dS > 5,]
outliers <- outliers[, c("abbrevs.x","mspec.x", "roi.x","roi.y", "dS")]
outliers %>% group_by(abbrevs.x,)
library(htmlTable)
library(magrittr)
outliers %>% htmlTable()



