# TO DO: bootstraps at end of dataset not yet integrated
## Series of 1way ANOVAS to test for lum diffs between the 3 subs each pop 
## Series of MANOVAS to test for col diffs betw substrates

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
df <- df %>% filter(substrate == "a") %>% group_by(abbrevs, mspec) %>% 
  sample_n(size = 1) %>% rbind(., df[df$substrate != "a",])

###* All dists 
jnds_c <- colDistEffic(bg_subs = "c") %>% 
  filter(abbrevs.x == abbrevs.y)
jnds_b <- colDistEffic(bg_subs = "b") %>%
  filter(abbrevs.x == abbrevs.y)

###* Nearest dists
# nearest 10% each substrate OR 
# nearest n ROIs across all images for a single plant in each image
nrstSub_dL_c <- nearSubset(jnds_c, localOnly = T, acrossImage = T, 
                           visfeature = "dL", qtiles =  10) 
nrstSub_dL_b <- nearSubset(jnds_b, localOnly = T, acrossImage = T, 
                           visfeature = "dL", qtiles =  10) 
nrstSub_dS_c <- nearSubset(jnds_c, localOnly = T, acrossImage = T, 
                           visfeature = "dS", qtiles =  10) 
nrstSub_dS_b <- nearSubset(jnds_b, localOnly = T, acrossImage = T, 
                           visfeature = "dS", qtiles =  10) 

homogDF <- function(df,dropInfo,substrate){ # format where substrate and visInfo can be disentangled 
  visInfoType <- c("dL", "dS")              # and kept in same df
  keepInfo <- visInfoType[!visInfoType %in% dropInfo]
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

# Mean distance metrics ####
###* Background-transfer metric ####
perLithMeans <- nrstSub_bc %>% 
  group_by(visInfo,subs,abbrevs.x, abbrevs.y, mspec.x, mspec.y) %>% 
  summarise(dAcrossImages = mean(distance)) %>%
  group_by(visInfo,subs,abbrevs.x, abbrevs.y, mspec.x) %>%
  summarise(dPerLith = mean(dAcrossImages))  # the mean of each image / lithops is one data point

perLithMeans <- perLithMeans %>% 
  mutate(collapseSpp.x = df$collapseSpp[match(abbrevs.x,df$abbrevs)], # import species & geol
         collapseSpp.y = df$collapseSpp[match(abbrevs.y,df$abbrevs)], 
         geology.x = df$geology[match(abbrevs.x,df$abbrevs)],
         geology.y = df$geology[match(abbrevs.y,df$abbrevs)]) %>%
  group_by(collapseSpp.x) %>% 
  mutate(popReps = n_distinct(abbrevs.x)) %>% 
  ungroup() %>%
  mutate(collapseSpp.x = fct_reorder(collapseSpp.x, popReps, min, .desc = T),
         collapseSpp.x = fct_relevel(collapseSpp.x, "singlePopSpp", after = Inf),
         visInfo = fct_recode(visInfo, LumContrast = "dL", ChromContrast = "dS")) %>%
  select(-popReps) %>% filter(abbrevs.x == abbrevs.y)

###** T-test & annotation df 
stat.test <- perLithMeans %>% mutate(subs = as_factor(subs)) %>%
  group_by(visInfo,abbrevs.x) %>% 
  t_test(dPerLith ~ subs) %>% # problem here, used t_test to make own f above
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() 
stat.test <- stat.test %>% add_xy_position(x = "abbrevs.x", scales = "free")
stat.test$collapseSpp.x <- perLithMeans$collapseSpp.x[match(stat.test$abbrevs.x,
                                                            perLithMeans$abbrevs.x)]
# boxplot

p <- ggplot(data = perLithMeans, 
            aes(x = abbrevs.x, y =  dPerLith)) +
  geom_boxplot(aes(fill = subs), outlier.shape = NA) + coord_flip() + 
  stat_pvalue_manual(stat.test, x = "abbrevs.x", coord.flip = T,
                     hide.ns = T, label = "p.adj.signif") +
  scale_fill_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                    labels = c("lithops-rock", "lithops-soil")) +
  theme_bw() + theme(panel.spacing = unit(0, "lines")) +
  labs(fill = "Substrate contrast", x = "Populations by species",
       y = "Noise scaled distance") +
  scale_y_continuous(sec.axis = dup_axis(name = waiver())) 
p <- p + facet_grid(vars(collapseSpp.x), vars(visInfo), scales = "free", space = "free_y") 
png(file = "output//rock-v-soil//boxplots-by-species-colour.png", type = 'cairo',
    width = 20, height = 30, units = "cm", res = 300, pointsize = 6) ; p; dev.off()

rm(stat.test, p, perLithMeans)
###* Paired distance metric ####
# Per image distances
source('functions//gmean.R')

# Nearest

nrstSub_bc <- nrstSub_bc %>% filter(mspec.x == mspec.y)

nrstLithops <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                   "abbrevs", "mspec","roi", "substrate")]), 
                     unique(nrstSub_bc[, c("abbrevs.x", "mspec.x", "roi.x")]),
                     by.y = c("abbrevs.x","mspec.x","roi.x"),
                     by.x = c("abbrevs","mspec", "roi")) %>%
  distinct(.) %>%
  select(-c('roi'))

nrstLumROIs <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                   "abbrevs", "mspec","roi", "substrate")]), 
                     unique(nrstSub_bc[
                       nrstSub_bc$visInfo == 'dL',
                       c("abbrevs.x", "mspec.x", "roi.y")]),
                     by.y = c("abbrevs.x","mspec.x","roi.y"),
                     by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(abbrevs, substrate, mspec) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean)

nrstChromROIs <- merge(unique(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                     "abbrevs", "mspec","roi", "substrate")]), 
                       unique(nrstSub_bc[
                         nrstSub_bc$visInfo == 'dS',
                         c("abbrevs.x", "mspec.x", "roi.y")]),
                       by.y = c("abbrevs.x","mspec.x","roi.y"),
                       by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(abbrevs, substrate, mspec) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean)

nrstLumLs <- rbind(nrstLithops, nrstLumROIs) %>%
  mutate(splitFct = paste(abbrevs,mspec, sep = '--')) %>%
  split(., .$splitFct)

nrstChromLs <- rbind(nrstLithops, nrstChromROIs) %>%
  mutate(splitFct = paste(abbrevs,mspec, sep = '--')) %>%
  split(., .$splitFct)

nrstLumLs <- lapply(nrstLumLs, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

nrstChromLs <- lapply(nrstChromLs, function(x) { # Pavo format
  rownames(x) <- make.unique(x$substrate)
  x <- rename(x, s = 'swMean', m = 'mwMean', l = 'lwMean', lum = 'lumMean')
  x <- x[,-which(!names(x) %in% c("s", "m", "l","lum"))]
})

nrstLumLs <- lapply(nrstLumLs, coldist,
                    n = c(1,16,32), weber = .05, weber.achro = .1, 
                    weber.ref = "longest", qcatch = "Qi", achromatic = TRUE)

nrstChromLs <- lapply(nrstChromLs, coldist, 
                      n = c(1,16,32), weber = .05, weber.achro = .1, 
                      weber.ref = "longest", qcatch = "Qi", achromatic = TRUE)

lumDists <- rbindlist(nrstLumLs, idcol= "names") %>%
  separate(., names, into = c("abbrevs", "mspec"), sep = "--") %>%
  select(-c(dS)) %>%
  filter(patch1 == 'a' | patch2 == 'a')

chromDists <- rbindlist(nrstChromLs, idcol= "names") %>%
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

rm(chromdis, lstLumDists, nrstChromROIs,
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
                           n = c(1,16,32), weber = .05, weber.achro = .1, 
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
               dPerLith = case_when(visInfo == 'dL' ~ 12.5,
                                    TRUE ~ 3.5),
               significanceDirection = case_when(estimate > 0 ~ 'Soil< Rock',
                                                 TRUE ~ 'Rock< Soil'),
               facetVar = paste(visInfo, dataset))

# Results 
# Turn this into tables
pairedTests %>%
  filter(dataset == 'Average', visInfo == 'dL', p.adj < 0.05, 
         significanceDirection == 'Rock< Soil') %>%
  count() /56

pairedTests %>%
  filter(dataset == 'Nearest 10%', visInfo == 'dL', p.adj < 0.05, 
         significanceDirection == 'Rock< Soil') %>%
  count() /56

pairedTests %>%
  filter(dataset == 'Average', visInfo == 'dL', p.adj < 0.05, 
         significanceDirection != 'Rock< Soil') %>%
  count() /56

pairedTests %>%
  filter(dataset == 'Nearest 10%', visInfo == 'dL', p.adj < 0.05, 
         significanceDirection != 'Rock< Soil') %>%
  count() /56

# mean JNDs < 1
inImageMeans %>% ungroup %>%
  filter(dataset == 'Average', visInfo == 'dL', c<1) %>%
  count()

###* boxlots
# boxplots idea stolen from https://www.stat.auckland.ac.nz/~paul/RGraphics/examples-dotplot.R

# nearest dist boxplots - inImageMeans
# all dist boxplots - inImageMeansAll
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
  geom_point(size = 0, stroke = 0)+
  scale_fill_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                    labels = c("lithops-rock", "lithops-soil"))+
  scale_colour_manual(values = c("Rock< Soil" = "#00BFC4",
                                 "Soil< Rock" = "#F8766D"))+
  guides(colour = guide_legend(
    override.aes = list(size = 5, shape = c(utf8ToInt("*"), utf8ToInt("*"))))) +
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

rm(pairedTests, p)

# LR V LS pop-mean distances (scatterplots) ####

source("functions\\gMeanDists.R")

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
                   n = c(1,16,32), weber = .05, weber.achro = .1, 
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
                     n = c(1,16,32), weber = .05, weber.achro = .1, 
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
                   n = c(1,16,32), weber = .05, weber.achro = .1, 
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


# Add substrate frequency variable
source("rock-v-soil//rock-estimate.R")
props <- df %>% count(abbrevs, substrate)
props <- props %>% pivot_wider(names_from = substrate, values_from = n)
props$rock_cover <- (props$b / (props$b + props$c)*100)
props$soil_cover <- (props$c / (props$b + props$c)*100)
props <- arrange(props, as.factor(abbrevs))
allGmeanDist <- merge(props, allGmeanDist, by.x = "abbrevs", by.y = "bgpop")
source("loading-cleaning.R") # not ideal. could change name of rock estimate df

allGmeanDist <- allGmeanDist %>% 
  mutate(collapseSpp = df$collapseSpp[match(.$abbrevs, df$abbrevs)],
         siteLabs = unlist(lapply(str_split(abbrevs, "_"), `[[`, 1)))

bindGmean$facetVar = paste(bindGmean$visfeature, bindGmean$pcentNrst, sep = '_')

bindGmean<- bindGmean %>% 
  mutate(pcentNrst = fct_recode(pcentNrst, '10%' = '10',
                                'All' = '100'))

p1<- ggplot(gMeans %>% filter(visInfo == "dL"), aes(x =c, y = b))+
  geom_point()+
  xlim(c(0,14))+ ylim(c(0,14))+
  xlab(NULL)+ylab('L-R gmean distance')+
  labs(title = "Luminance")+
  geom_abline(lty = 2)+
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
  caption = 'L-S gmean distance',
  theme = theme(
    plot.caption = element_text(hjust = .55,margin = margin(t=0), size = 11.5))) +
  plot_layout(ncol = 2)
dev.off()

rm(p,p1,p2)

# Results
bindGmean %>%
  filter(pcentNrst =='10%', visfeature == 'dL', b<1 | c<1 ) 
35/56
bindGmean %>%
  filter(pcentNrst =='All', visfeature == 'dL', b<1 | c<1 ) 
13/56

bindGmean %>%
  filter(pcentNrst =='All', visfeature == 'dL', b<1) %>%
  count()/56

bindGmean %>%
  filter(pcentNrst =='All', visfeature == 'dL', b<1 | c<1 ) 


bindGmean %>%
  filter(pcentNrst =='10%', visfeature == 'dS', b<1) 
bindGmean %>%
  filter(pcentNrst =='All', visfeature == 'dS', b<1) %>%
  count()

bindGmean %>%
  filter(pcentNrst =='10%', visfeature == 'dS', c<1) %>%
  count()
bindGmean %>%
  filter(pcentNrst =='All', visfeature == 'dS', c<1) %>%
  count()

# 
# ###* Scatterplots LR vs LS (average dataset)
# dummy <- data.frame('a_b' =  c(range(allGmeanDist$dL_a_b,
#                                       allGmeanDist$dL_a_c),
#                                range(allGmeanDist$dS_a_b,
#                                      allGmeanDist$dS_a_c)),
#                     'a_c' = c(range(allGmeanDist$dL_a_b,
#                                   allGmeanDist$dL_a_c),
#                               range(allGmeanDist$dS_a_b,
#                                     allGmeanDist$dS_a_c)),
#                     'visInfo' = c(rep('dL',2), rep('dS',2)),
#                     'rock_cover' = 0
#                     )
# png("output//rock-v-soil//LR-v-LS-facet-gMeanDist-scatter.png",
#     type = 'cairo', units = 'mm', width = 215, height = 300, res = 300)
# allGmeanDist %>% pivot_longer(cols = contains(c('dL', 'dS')),
#                               names_to = c('visInfo','.value'),
#                               names_sep = 3) %>%
#   mutate(visInfo = fct_recode(visInfo, 'dL' = 'dL_', 'dS' = 'dS_')) %>%
#   ggplot(aes(x = a_c, y = a_b, colour = rock_cover))+ 
#   geom_point()+
#   geom_abline(slope = 1, intercept = 0, lty = 2)+ 
#   geom_blank(data = dummy)+
#   theme_bw()+
#   #coord_fixed(xlim = c(0,2.2), ylim = c(0,2.2))+ 
#   labs(y = "L-R geometric mean distance (JNDs)",
#        x = "L-S geometric mean distance (JNDs)")+
#   theme(axis.title = element_text(size = 10), 
#         aspect.ratio = 1)+
#   facet_wrap(vars(visInfo), scales = 'free', ncol = 1,
#              labeller = labeller(visInfo = c(dL = 'Lum dist', dS = 'Chroma dist')))
# dev.off()
#   
# 
# ###** small multiples facet by spp
# pLum <- ggplot(allGmeanDist, aes(x = dL_a_c, y = dL_a_b))+
#   geom_point(data = allGmeanDist[, c("dL_a_c", "dL_a_b")], color = "gray")+
#   geom_point(aes(fill = rock_cover), size = 2, pch =21)+ 
#   geom_abline(intercept = 0, lty = 2)+
#   geom_text_repel(aes(label = siteLabs),
#                   max.overlaps = Inf,
#                   min.segment.length = 0,
#                   box.padding = 0.5)+
#   xlab("LS geometric mean distance")+ ylab("LR geometric mean distance")+
#   facet_wrap(vars(collapseSpp))+
#   theme(aspect.ratio = 1)
# png("output\\rock-v-soil\\LR-vs-LS-lum-gmeandist-scatter_facetSpp.png",
#     type = 'cairo',
#     width = 150, height = 150, units = "mm", res = 300, pointsize = 6)
# pLum; dev.off()
# 
# pCol <- ggplot(allGmeanDist, aes(x = dS_a_c, y = dS_a_b))+
#   geom_point(data = allGmeanDist[, c("dS_a_c", "dS_a_b")], color = "gray")+
#   geom_point(aes(fill = rock_cover), size = 2, pch =21)+ 
#   geom_abline(intercept = 0, lty = 2)+
#   geom_text_repel(aes(label = siteLabs),
#                   max.overlaps = Inf,
#                   min.segment.length = 0,
#                   box.padding = 0.5)+
#   xlab("LS geometric mean distance")+ ylab("LR geometric mean distance")+
#   facet_wrap(vars(collapseSpp))+
#   theme(aspect.ratio = 1)
# png("output\\rock-v-soil\\LR-vs-LS-col-gmeandist-scatter_facetSpp.png",
#     type = 'cairo',
#     width = 150, height = 150, units = "mm", res = 300, pointsize = 6)
# pCol; dev.off()
# 
# rm(pCol, pLum)

###** Overall Boxplots pop-mean LR , LS dists ####
### paired t-test difference in overall LR-LS dists

overallTtest <- bindGmean %>%
  group_by(pcentNrst, visfeature) %>%
  nest() %>%
  mutate(ttest = map(data, t_test)) %>%
  unnest(ttest) %>% ungroup()

overallTtest <- overallTtest %>%
  mutate(labelz = 
           case_when(p.value < 0.001 ~ '***',
                     p.value < 0.01 ~ '**',
                     p.value < 0.05 ~ '*'),
         meanDist = case_when(visfeature == 'dL' ~ 9,
                          visfeature == 'dS' ~ 1.75))

pData <- bindGmean %>%
  pivot_longer(cols = c("b","c"),
               names_to = "substrate", values_to = "meanDist")
p2 <-ggplot(pData, aes(pcentNrst,meanDist))+
  geom_boxplot(aes(fill = substrate))+
  geom_text(data = overallTtest, aes(label = labelz), size = 5)+
  theme_classic()+
  labs(y = 'Population-mean distances',
       fill = 'Substrate contrast')+
  scale_x_discrete('Dataset', labels = c('10%' = 'Nearest 10%',
                                         'All' = 'Average') )+
  scale_fill_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                    labels = c("lithops-rock", "lithops-soil"))+
  facet_grid(vars(visfeature), scales = "free_y",
             labeller = labeller(visfeature = c('dL' = 'Luminance',
                                                'dS' = 'Chroma')))


png("output\\rock-v-soil\\LR-LS-dist-relship-varyNrstPcent-boxplots.png",
    type = 'cairo',units = 'mm', res = 300, width = 200, height = 140)
p2
dev.off()

# Results


rm(pdata, overallTtest)

###** STRANGE RESULT ####
# TK verruc HG smaller distance L-S for the all subset because it has a 
# different subset of Lithops. Because there are only a few images with soils, 
# the nearest subsets calculate the location of only the few Lithops which are
# in those images  

###** Table summary stats for nearest sets of dists ####
nMspecs <- df %>%
  group_by(abbrevs) %>%
  summarise(nMspec = n_distinct(mspec))

tenPcent <- rbindlist(distSubsets, idcol = T) %>%
  filter((.id == "c_10%_dL" | .id == "b_10%_dL") & substrate != "a")
maxCount <- tenPcent %>%
  group_by(abbrevs, substrate, mspec) %>%
  summarise(nMax = max(n())) %>%
  group_by(abbrevs, substrate) %>%
  summarise(nMax = max(nMax)) %>% 
  pivot_wider(names_from = substrate, values_from = nMax, 
              names_glue  = "{substrate}_{.value}")

tenPcent <- tenPcent %>% 
  group_by(abbrevs, substrate) %>%
  summarise(nAll = n(),
            noImages = n_distinct(mspec)) %>%
  pivot_wider(names_from = substrate, values_from = c(noImages, nAll), 
              names_glue  = "{substrate}_{.value}")

tenPcent <- merge(tenPcent, nMspecs, by = "abbrevs") %>%
  merge(.,maxCount) 
tenPcent <- tenPcent %>% relocate(nMspec,.after = abbrevs)

rm(tenPcent, nMspecs, maxCount)
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
# To do - lumMean isn't the perceptually transformed version, nor is it log scaled etc.


###* plotting dfs ####

nrstLithops <- merge(unique(df[, c('xCoord', 'yCoord',"lumMean",
                                   "abbrevs", "mspec","roi", "substrate")]), 
                     unique(nrstSub_bc[, c("abbrevs.x", "mspec.x", "roi.x")]),
                     by.y = c("abbrevs.x","mspec.x","roi.x"),
                     by.x = c("abbrevs","mspec", "roi")) %>%
  distinct(.) %>%
  select(c(abbrevs,substrate, xCoord, yCoord, lumMean))

gMeanLum <- merge(unique(df[, c('xCoord', 'yCoord',"lumMean",
                                "abbrevs", "mspec","roi", "substrate")]), 
                  unique(nrstSub_bc[
                    nrstSub_bc$visInfo == 'dL',
                    c("abbrevs.x", "mspec.x", "roi.y")]),
                  by.y = c("abbrevs.x","mspec.x","roi.y"),
                  by.x = c("abbrevs","mspec", "roi")) %>%
  group_by(substrate,abbrevs, mspec) %>% 
  summarise_at(c("xCoord", "yCoord", "lumMean"), mean) %>%
  group_by(substrate, abbrevs) %>% # this step not all that NB, highly correlated w straight mean
  summarise_at(c("xCoord", "yCoord", "lumMean"), mean) %>%
  mutate
  pivot_wider(names_from = substrate,
              values_from = c(xCoord, yCoord, lumMean))
  
gMeanChrom <- merge(unique(df[, c("xCoord", "yCoord","lumMean",
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

distSubsetsTf <- distSubsets %>% 
  rbindlist(idcol = T) %>%
  merge(., df[, c(colnames(.)[colnames(.) %in% colnames(df)],
                     "xCoord", "yCoord")]) %>%
  group_by(.id,substrate,abbrevs, mspec) %>% 
  summarise_at(c("xCoord", "yCoord", "lumMean"), mean) %>%
  group_by(.id,substrate, abbrevs) %>% # this step not all that NB, highly correlated w straight mean
  summarise_at(c("xCoord", "yCoord", "lumMean"), mean) %>%
  separate(., col = .id,
           into = c("substratee","nrstSet", "visInf")) %>%
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
  mutate(nrstSet = as_factor(nrstSet),
         nrstSet = fct_relevel(nrstSet, "100", after = Inf))

x <- distSubsetsTf %>%
  pivot_longer(cols = c(contains(c("Rock","Soil")))) %>% 
  separate(., "name",into = c("predVisInfo", "predSub"), sep = "_")
###* Scatterplot L coords as function R or S coords (all points) ####
library(ggpmisc)
library(ggrepel)

meanCoords <- x %>% pivot_longer(cols = contains(c("lithops")),
                             names_to = c("lithVisInfo", "lithSub"),
                             names_sep = "_",
                             values_to = 'lithValue') %>%
  mutate(facetVar = paste(predVisInfo,predSub)) %>%
  filter(predVisInfo == lithVisInfo, nrstSet == "100") %>%
  select(-c(lithSub, lithVisInfo)) %>%
  cbind(.,props[match(.$abbrevs, props$abbrevs),
                         c("rock_cover", "soil_cover")]) %>%
  cbind(., df[match(.$abbrevs, df$abbrevs),
              c("collapseSpp", "pop_spp")]) %>%
  separate(., abbrevs, "siteLabs", extra = 'drop')

#dummy <- data.frame(range)

facetVarRename = c('xCoord Rock' = 'GR Rock', 'xCoord Soil' = 'GR Soil',
                   'yCoord Rock' = 'BY Rock', 'yCoord Soil' = 'BY Soil',
                   'lumMean Rock' = 'Lum Rock', 'lumMean Soil' = 'Lum Soil')
p <- ggplot(meanCoords,aes(value,lithValue, colour = rock_cover))+
  geom_point()+
  geom_abline(lty='dashed')+
  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'black', se = F, lwd = 0.1)+
  facet_wrap(vars(facetVar), scales = 'free',
             ncol = 2, strip.position = 'bottom',
             labeller = labeller(facetVar = facetVarRename))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        aspect.ratio = 1,
        strip.placement = 'outside',
        strip.background = element_rect(fill = NA, color = NA))+
  xlab(NULL)+
  ylab("Lithops mean")+
  stat_poly_eq(formula = y~x,
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                 sep = "~~~")), parse = TRUE, vstep = 3)
 
png("output\\rock-v-soil\\L-coords-af-RS-coords-scatter.png", type = 'cairo',
    units = "mm", res = 300, width = 215, height = 300)
p
dev.off()

###** Small multiples facet by spp ####

meanCoords <- x %>% pivot_longer(cols = contains(c("lithops")),
                                 names_to = c("lithVisInfo", "lithSub"),
                                 names_sep = "_",
                                 values_to = 'lithValue') %>%
  mutate(facetVar = paste(predVisInfo,predSub)) %>%
  filter(predVisInfo == lithVisInfo, nrstSet == "100") %>%
  select(-c(lithSub, lithVisInfo)) %>%
  cbind(.,props[match(.$abbrevs, props$abbrevs),
                c("rock_cover", "soil_cover")]) %>%
  cbind(., df[match(.$abbrevs, df$abbrevs),
              c("collapseSpp", "pop_spp")]) %>%
  separate(., abbrevs, "siteLabs", extra = 'drop')

library(ggforce)
library(cowplot)

sppLabeller<- c(olivacea = 'L. olivacea', localis = 'L. localis',
                marmorataSpp = 'L. marmorata', divergensSpp = 'L. divergens',
                halliiSpp = 'L. hallii', leslieiSpp = 'L. lesliei',
                otzeniana = 'L. otzeniana', compt. = 'L. comptonii',
                fulleri = 'L. fulleri', bromf. = 'L. bromfieldii',
                hookeriSpp = 'L. hookeri', meyeri = 'L. meyeri', 
                dint. = 'L. dinteri', verruc. = 'L. verruculosa', 
                singlePops = 'singlePops')
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
p2 <-vector(mode = 'list', length = 2)
p3 <- vector(mode = 'list', length = 2)

for(i in 1:2){
  j <- seq(1, length(unique(meanCoords$collapseSpp)), 10)
  
  pDataBg<-meanCoords %>%
    filter(nrstSet == '100',
           predVisInfo == 'lumMean')
  pDataOver <- meanCoords[
    meanCoords$collapseSpp %in% 
      unique(meanCoords$collapseSpp)[j[i]:(j[i]+9)],] %>%
    filter(nrstSet == '100',
           predVisInfo == 'lumMean')
  
  p1[[i]]<- ggplot(pDataOver,aes(x= value, y=lithValue))+
    geom_point(data = pDataBg[, c("predSub","value", "lithValue")],
               color = "gray")+
    geom_point(aes(fill = rock_cover), size = 2, pch =21)+ 
    geom_abline(intercept = 0, lty = 2)+
    #labs(title = 'Luminance')+
    xlab('Substrate luminance')+
    ylab('Lithops luminance')+
    facet_grid(collapseSpp~predSub, switch = 'y', 
               labeller = labeller(collapseSpp = sppLabeller))+
    custText+
    theme(strip.placement = 'outside',
          strip.text.y = element_text(face = 'italic'),
          strip.background = element_blank(),
          #panel.spacing.y = unit(0,"lines"),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none')
  
  pDataBg<-meanCoords %>%
    filter(nrstSet == '100',
           predVisInfo == 'xCoord')
  pDataOver <- meanCoords[
    meanCoords$collapseSpp %in% 
      unique(meanCoords$collapseSpp)[j[i]:(j[i]+9)],] %>%
    filter(nrstSet == '100',
           predVisInfo == 'xCoord')
  
  p2[[i]]<- ggplot(pDataOver,aes(x= value, y=lithValue))+
    geom_point(data = pDataBg[, c('predSub',"value", "lithValue")],
               color = "gray")+
    geom_point(aes(fill = rock_cover), size = 2, pch =21)+ 
    geom_abline(intercept = 0, lty = 2)+
    custText+
    #labs(title = 'G:R opponent')+
    xlab('Substrate G:R opponent')+
    ylab('Lithops G:R opponent')+
    facet_grid(collapseSpp~predSub)+
    custTheme
  
  pDataBg<-meanCoords %>%
    filter(nrstSet == '100',
           predVisInfo == 'yCoord')
  pDataOver <- meanCoords[
    meanCoords$collapseSpp %in% 
      unique(meanCoords$collapseSpp)[j[i]:(j[i]+9)],] %>%
    filter(nrstSet == '100',
           predVisInfo == 'yCoord')
  
  p3[[i]]<- ggplot(pDataOver,aes(x= value, y=lithValue))+
    geom_point(data = pDataBg[, c('predSub',"value", "lithValue")],
               color = "gray")+
    geom_point(aes(fill = rock_cover), size = 2, pch =21)+ 
    geom_abline(intercept = 0, lty = 2)+
    custText+
    #labs(title = 'Y:B opponent')+
    xlab('Substrate Y:B opponent')+
    ylab('Lithops Y:B opponent')+
    facet_grid(collapseSpp~predSub)+
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

# smallMultPlots <- vector(mode = 'list', length =  6)
# names(smallMultPlots) <- c('pDataLumR', 'pDataLumS','pDataxCoordR',
#                            'pDataxCoordS','pDatayCoordR', 'pDatayCoordS')
# smallMultPlots[[1]] <- meanCoords %>% filter(predVisInfo == 'lumMean' & 
#                                  predSub == 'Rock')
# smallMultPlots[[2]] <- meanCoords %>% filter(predVisInfo == 'lumMean' & 
#                                      predSub == 'Soil')
# smallMultPlots[[3]]  <- meanCoords %>% filter(predVisInfo == 'xCoord' & 
#                                      predSub == 'Rock')
# smallMultPlots[[4]]  <- meanCoords %>% filter(predVisInfo == 'xCoord' & 
#                                         predSub == 'Soil')
# smallMultPlots[[5]]  <- meanCoords %>% filter(predVisInfo == 'yCoord' & 
#                                         predSub == 'Rock')
# smallMultPlots[[6]]  <- meanCoords %>% filter(predVisInfo == 'yCoord' & 
#                                         predSub == 'Soil')
# 
# plotF <- function(pData){
#   
#   ggplot(pData, aes(value, lithValue))+
# geom_point(data = pData[, c("value", "lithValue")],
#            color = "gray")+
#   geom_point(aes(fill = rock_cover), size = 2, pch =21)+ 
#   geom_abline(intercept = 0, lty = 2)+
#   geom_text_repel(aes(label = siteLabs),
#                   force = .5,
#                   size = 3,
#                   max.overlaps = Inf,
#                   min.segment.length = 0.5)+
#   xlab(unique(paste(pData$predVisInfo, pData$predSub)))+
#   ylab(unique(paste(pData$predVisInfo, "Lithops")))+
#   facet_wrap(vars(collapseSpp), ncol = 5)+
#   theme(aspect.ratio = 1, legend.position = 'none')
# }

smallMultPlots <- lapply(smallMultPlots, plotF)
lapply(names(smallMultPlots), 
       function(x) 
       ggsave(paste0("output\\rock-v-soil\\L-coords-af-RS-coords-scatter",
       x,".png"), plot=smallMultPlots[[x]],
       units = 'mm', width = 215, height = 140, dpi = 300))

rm(smallMultPlots)
###* Scatterplot L coords as function R or S coords (only dominant) ####

png("output\\rock-v-soil\\L-coords-af-RS-coords-scatter-only-dominant.png",
    type = 'cairo', units = "mm", res = 300, width = 215, height = 300)
meanCoords %>% filter((predSub == 'Rock' & rock_cover > 70) | 
                        predSub == 'Soil' & soil_cover > 60) %>%
ggplot(aes(value,lithValue, colour = rock_cover))+
  geom_point()+
  geom_abline(lty='dashed')+
  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'black', se = F, lwd = 0.1)+
  facet_wrap(vars(facetVar), scales = 'free',
             ncol = 2, strip.position = 'bottom')+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        aspect.ratio = 1,
        strip.placement = 'outside',
        strip.background = element_rect(fill = NA, color = NA))+
  xlab(NULL)+
  ylab("Lithops mean")+
  stat_poly_eq(formula = y~x,
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                 sep = "~~~")), parse = TRUE, vstep = 3)
dev.off()

###* Scatterplot L coords as function R or S coords varying nearest % ####

meanCoords <- x %>% pivot_longer(cols = contains(c("lithops")),
                                 names_to = c("lithVisInfo", "lithSub"),
                                 names_sep = "_",
                                 values_to = 'lithValue') %>%
  filter(predVisInfo == lithVisInfo, nrstSet != '50') %>%
  select(-c(lithSub, lithVisInfo)) %>%
  cbind(., df[match(.$abbrevs, df$abbrevs),
              c("collapseSpp", "pop_spp")]) %>%
  separate(., abbrevs, "siteLabs", extra = 'drop')

nrstSetLabeller <- c('10' = 'Nearest 10%', '100' = 'Average')
predVisInfoLabeller <- c('lumMean' = 'Luminance',
                         xCoord = 'G:R opponent', yCoord = 'Y:B opponent')
statPoly <- stat_poly_eq(
  formula = y~x, 
  aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label..,
                    sep = "~~~")), parse = TRUE, vstep = 3, size =3)
# custTheme <- theme(strip.background = element_blank(),
#       strip.placement = 'outside',
#       strip.text.x = element_text(margin = margin(2, 0, 2, 0)),
#       strip.text.y = element_blank(),
#       ggh4x.facet.nestline = element_line(),
#       axis.title.y = element_text(size = rel(.8), angle = 90))
# custFacet <- facet_nested(predVisInfo ~ predSub +nrstSet,
#                         #switch = 'x',
#                         nest_line = element_line(),
#                         labeller = labeller(nrstSet = nrstSetLabeller,
#                                             predVisInfo = predVisInfoLabeller))

library(patchwork)

p1<-meanCoords %>% 
  filter(predVisInfo == 'lumMean') %>%
  ggplot(aes(x = value, y = lithValue))+
  geom_point()+
  geom_abline(lty='dashed')+
  geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .7)+
  #coord_equal()+
  #xlab(NULL)+
  xlab('Background subs. luminance')+
  ylab('Lithops luminance')+
  statPoly+ 
  facet_grid(nrstSet~predSub,
             labeller = labeller(nrstSet = nrstSetLabeller))+
  theme(strip.background = element_rect(fill = 'grey91'))

p2<-meanCoords %>% 
  filter(predVisInfo == 'xCoord') %>%
  ggplot(aes(x = value, y = lithValue))+
  geom_point()+
  geom_abline(lty='dashed')+
  geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .7)+
  #coord_equal()+
  #xlab(NULL)+
  xlab('Background subs. G:R opponent')+
  ylab('Lithops G:R opponent')+
  statPoly+ 
  facet_grid(nrstSet~predSub,
             labeller = labeller(nrstSet = nrstSetLabeller))+
  theme(strip.background = element_rect(fill = 'grey91'))

p3<-meanCoords %>% 
  filter(predVisInfo == 'yCoord') %>%
  ggplot(aes(x = value, y = lithValue))+
  geom_point()+
  geom_abline(lty='dashed')+
  geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .7)+
  #coord_equal()+
  #xlab(NULL)+
  xlab('Background subs. Y:B opponent')+
  ylab('Lithops Y:B opponent')+
  statPoly+ 
  facet_grid(nrstSet~predSub,
             labeller = labeller(nrstSet = nrstSetLabeller))+
  theme(strip.background = element_rect(fill = 'grey91'))


png("output\\rock-v-soil\\L-coords-af-RS-coords-amalg.png",
    type = 'cairo', units = "mm", res = 300,
    width = 215, height = 260)
p1+p2+p3 + plot_layout(nrow = 3)
dev.off()
rm(custFacet, custTheme,statPoly)

# p1<-meanCoords %>% 
#   filter(predVisInfo == 'lumMean') %>%
# ggplot(aes(x = value, y = lithValue))+
#   geom_point()+
#   geom_abline(lty='dashed')+
#   geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .5)+
#   coord_equal()+
#   #xlab(NULL)+
#   xlab('Substrate luminance')+
#   ylab('Lithops luminance')+
#   statPoly+custTheme+custFacet
# 
# p2<- meanCoords %>% 
#   filter(predVisInfo == 'xCoord') %>%
#   ggplot(aes(x = value, y = lithValue))+
#   geom_point()+
#   geom_abline(lty='dashed')+
#   geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .5)+
#   coord_equal()+
#   #xlab(NULL)+
#   xlab('Substrate G:R opponent')+
#   ylab('Lithops G:R opponent')+
#   statPoly+custTheme+custFacet
# 
# p3<- meanCoords %>% 
#   filter(predVisInfo == 'yCoord') %>%
#   ggplot(aes(x = value, y = lithValue))+
#   geom_point()+
#   geom_abline(lty='dashed')+
#   geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .5)+
#   coord_equal()+
#   #xlab(NULL)+
#   xlab('Substrate Y:B opponent')+
#   ylab('Lithops Y:B opponent')+
#   statPoly+custTheme+custFacet



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
                                          weber = .1,
                                          weber.achro = .1, 
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

# Rock spread predicts Lithops distance to rock ####
source("functions//MAD.R")
distSpreadAll <- df %>% group_by(abbrevs, substrate) %>% 
  summarize(MADXcoord = MAD(xCoord),
            MADYcoord = MAD(yCoord),
            MADLum = MAD(lumMean)) %>% 
  pivot_wider(names_from = "substrate",
              values_from = c("MADXcoord", "MADYcoord", 'MADLum')) %>%
  left_join(allGmeanDist, .)

library(patchwork)

p1 <- distSpreadAll %>%
  filter(rock_cover > 70) %>%
  ggplot(aes(MADLum_b, dL_a_b))+
  geom_point()+
  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'black', se = F, lwd = 0.5)+
  ylab('LR luminance gMean dist')+
  xlab('Luminance MAD')+
  stat_poly_eq(formula = y~x,
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                 sep = "~~~")), parse = TRUE, vstep = 3,
               size = 2.5)+
  theme_classic()

p2 <- distSpreadAll %>%
  filter(rock_cover > 70) %>%
  ggplot(aes(MADXcoord_b, dS_a_b))+
  geom_point()+
  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'black', se = F, lwd = 0.5)+
  ylab('LR chromatic gMean dist')+
  xlab('RG MAD')+
  stat_poly_eq(formula = y~x,
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                 sep = "~~~")), parse = TRUE, vstep = 3,
               size =2.5)+
  theme_classic()

p3 <- distSpreadAll %>%
  filter(rock_cover > 70) %>%
  ggplot(aes(MADYcoord_b, dS_a_b))+
  geom_point()+
  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'black', se = F, lwd = 0.5)+
  ylab('LR chromatic gMean dist')+
  xlab('BY MAD')+
  stat_poly_eq(formula = y~x,
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                 sep = "~~~")), parse = TRUE, vstep = 3,
               size = 2.5)+
  theme_classic()
png(file = "output//rock-v-soil//R-spread-pred-LR-dist.png", type = 'cairo',
    width = 215, height = 140, units = "mm", res = 300) 
p1 + p2 + p3 + plot_layout(nrow = 1)
dev.off()

rm(distSpreadAll, p1,p2,p3)

# Other stuff etc. ####

# outliers
outliers <- jnds_dec[jnds_dec$dS > 5,]
outliers <- outliers[, c("abbrevs.x","mspec.x", "roi.x","roi.y", "dS")]
outliers %>% group_by(abbrevs.x,)
library(htmlTable)
library(magrittr)
outliers %>% htmlTable()

# 2: stat tests
## paired t-tests on L-R and L-S distances for substrate matching
### Luminance - nearest 10% 
library(broom)
library(purrr)
batchTest <- function(df){with(df, by(df, abbrevs.x, 
                                      function(x) t.test(meanDistPerLith ~ bgSubstrate, data=x)))}

nrstJnds <- nrstSub # take dataset computed above

perLithMeans <- nrstJnds %>% 
  group_by(abbrevs.x, bgSubstrate, mspec.x, mspec.y) %>% 
  summarise(dLAcrossImages = mean(dL)) %>%
  group_by(abbrevs.x, bgSubstrate, mspec.x) %>%
  summarise(meanDistPerLith = mean(dLAcrossImages)) # the mean of each image / lithops is one data point

dLres <- batchTest(perLithMeans)
dLres <- map_df(dLres, tidy, .id = "abbrevs.x")
dLres$species <- df$species[match(dLres$abbrevs.x,df$abbrevs)]
dLres <- dLres %>% 
  mutate(pAdjusted = p.adjust(dLres$p.value)) %>%
  mutate_if(is.numeric, signif, digits=2) %>%
  arrange(species) %>%
  select(-one_of("species", "estimate", "method", "alternative",
                 "conf.low", "conf.high")) %>%
  rename("Population" = "abbrevs.x", "LR mean dist" = "estimate1",
         "LS mean dist" = "estimate2", "t" = "statistic", "p" = "p.value",
         "df" = "parameter")
write.csv(dLres, quote = F)

### Colour
jndsbc <- rbind(jnds_b,jnds_c) # calculated above

perLithMeans <- jndsbc %>% 
  filter(abbrevs.x == abbrevs.y) %>%
  mutate(bgSubstrate = substr(roi.y, 0,1)) %>%
  group_by(abbrevs.x, bgSubstrate, mspec.x, mspec.y) %>% 
  summarise(dSAcrossImages = mean(dS)) %>%
  group_by(abbrevs.x, bgSubstrate, mspec.x) %>%
  summarise(meanDistPerLith = mean(dSAcrossImages)) # the mean of each image / lithops is one data point

dSres <- batchTest(perLithMeans)
dSres <- map_df(dSres, tidy, .id = "abbrevs.x")
dSres$species <- df$species[match(dSres$abbrevs.x,df$abbrevs)]
dSres <- dSres %>%
  mutate(pAdjusted = p.adjust(dSres$p.value)) %>%
  mutate_if(is.numeric, signif, digits=2) %>%
  arrange(species) %>% 
  select(-one_of("species", "estimate", "method", "alternative",
                 "conf.low", "conf.high")) %>%
  rename("Population" = "abbrevs.x", "LR mean dist" = "estimate1",
         "LS mean dist" = "estimate2", "t" = "statistic", "p" = "p.value",
         "df" = "parameter") 

write.csv(dSres, quote = F)


