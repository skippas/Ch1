# TO DO: bootstraps at end of dataset not yet integrated
## Series of 1way ANOVAS to test for lum diffs between the 3 subs each pop 
## Series of MANOVAS to test for col diffs betw substrates

# Matching rock or soil? T-tests & boxplots ####
# SETUP ####
library(ggpubr)
library(rstatix)
source("loading-cleaning.R")
source("functions\\nearSubset.R") 
source("functions\\coldist_effic.R")
source("functions\\distToOrigDF.R")

theme_set(theme_bw())
theme_update(panel.spacing = unit(0, "lines"),
             strip.placement = "outside",
             strip.background = element_rect(colour = "black")) # see riffomonas proj


# Generate distances ####
# A single Lithops and local comps only 
oneLith <- df %>% filter(substrate == "a") %>% group_by(abbrevs, mspec) %>% 
  sample_n(size = 1) 
df <- rbind(oneLith, df[df$substrate != "a",])
rm(oneLith)

# nearest 10% each substrate OR nearest n ROIs across all images for a single plant 
# in each image
jnds_c <- colDistEffic(bg_subs = "c")
jnds_b <- colDistEffic(bg_subs = "b")
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

###** T-test & annotation df ####
stat.test <- perLithMeans %>% mutate(subs = as_factor(subs)) %>%
  group_by(visInfo,abbrevs.x) %>% 
  t_test(dPerLith ~ subs) %>% # problem here, used t_test to make own f above
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() 
stat.test <- stat.test %>% add_xy_position(x = "abbrevs.x", scales = "free")
# This code below is a work around for not using it as a grouping var above
stat.test$collapseSpp.x <- perLithMeans$collapseSpp.x[match(stat.test$abbrevs.x,
                                                            perLithMeans$abbrevs.x)]
#stat.test <- stat.test %>% mutate(x.position = case_when(visInfo == "dS" ~ 2.5,
#                     visInfo == "dL" ~ 15))



###* Paired distance metric  ####
inImageMeans <- nrstSub_bc %>% filter(mspec.x == mspec.y) %>%
  group_by(visInfo,subs,abbrevs.x, abbrevs.y, mspec.x, mspec.y) %>% 
  summarise(dInImages = mean(distance)) %>% 
  pivot_wider(names_from = subs, values_from = dInImages)
inImageMeans$collapseSpp.x <- 
  df$collapseSpp[match(inImageMeans$abbrevs.x, df$abbrevs)]

###** Paired T-test & annotation df ####
t_test <- function(df, mu = 0, alt = "two.sided", paired = T,
                   conf.level = .95,var.equal = F){
  tidy(t.test(df$b, df$c,
              mu = mu, alt = alt,
              conf.level = conf.level,
              paired = paired, var.equal = var.equal))
}

pairedTests <- inImageMeans %>%
  group_by(abbrevs.x, visInfo) %>%
  nest() %>%
  mutate(ttest = map(data, t_test)) %>%
  unnest(ttest) %>% ungroup()

pairedTests[pairedTests$visInfo == "dL","p.adj"] <- 
  p.adjust(pairedTests[pairedTests$visInfo == "dL",]$p.value,
           n = 56, method = "bonferroni")
pairedTests[pairedTests$visInfo == "dS","p.adj"] <- 
  p.adjust(pairedTests[pairedTests$visInfo == "dS",]$p.value,
           n = 56, method = "bonferroni")
pairedTests <- pairedTests %>% mutate(labelz = case_when(p.value < 0.05 ~ "*"))
pairedTests$collapseSpp.x <- df$collapseSpp[match(pairedTests$abbrevs.x,df$abbrevs)]
pairedTests[pairedTests$visInfo == "dL", "dPerLith"] <- 12.5
pairedTests[pairedTests$visInfo == "dS", "dPerLith"] <- 3.3

###* boxlots #### 
# boxplots idea stolen from https://www.stat.auckland.ac.nz/~paul/RGraphics/examples-dotplot.R

# paired, nearest 10%? 
p <- inImageMeans %>%
  pivot_longer(cols = c(b,c), 
               names_to = "subs",
               values_to = "dPerLith") %>%
  ggplot(aes(x = abbrevs.x, y =  dPerLith))+
  geom_boxplot(aes(fill = subs))+
  geom_text(data = pairedTests,
            aes(label = labelz, x = abbrevs.x, y = dPerLith))+
  scale_fill_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                    labels = c("lithops-rock", "lithops-soil"))+
  labs(fill = "Substrate contrast",
       x = "Populations by species",
       y = "Noise scaled distance")+
  scale_y_continuous(sec.axis = dup_axis(name = waiver()))+
  coord_flip()+ 
  facet_grid(vars(collapseSpp.x), vars(visInfo),switch = "y",
             scales = "free", space = "free_y") 
png(file = "output//rock-v-soil//RvS-pairedTest.png", type = 'cairo',
    width = 20, height = 30, units = "cm", res = 300, pointsize = 6) ; p; dev.off()

rm(inImageMeans, pairedTests)
# background transfer

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

rm(stat.test, perLithMeans)

# ROCK V SOIL mean distance (scatterplots) ####
#* Plotting data ####
## luminance (the mean of the coordinates of the nearest set of distances)
source("functions\\gMeanDists.R")
nearJnds <- nrstSub_bc[nrstSub_bc$visInfo != "dS",]
nrDistDF <- distToOrigDF(distDF = nearJnds, origDF = df) 

nrDistDF <- nrDistDF %>% 
  group_by(abbrevs, substrate, mspec) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), mean) # each image contributes equally to gmean dist
nrGmeanDist <- gMeanDists(df = nrDistDF, combine_bg = F,substrates = c("b","c")) # gmeans

nrGmeanDist <- nrGmeanDist %>% 
  filter(comparison == "local") %>% 
  pivot_wider(values_from = c("dL","dS"),
              names_from = c("patch1", "patch2"))
rm(nearJnds, nrDistDF)

## Colour
df2 <- df %>% 
  group_by(abbrevs, substrate, mspec) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), mean) # each image contributes equally to gmean dist

allGmeanDist <- gMeanDists(df = df2, combine_bg = F,substrates = c("b","c")) %>%
  filter(comparison == "local") %>% 
  pivot_wider(values_from = c("dL","dS"),
              names_from = c("patch1", "patch2"))

#* Scatterplots ####
p <- ggplot(nrGmeanDist, aes(x = dL_a_c, y = dL_a_b)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + coord_fixed(xlim = c(0,15), ylim = c(0,15)) + 
  labs(y = "L-R geometric mean lum dist (nearest subset) (JNDs)",
       x = "L-S geometric mean lum dist (nearest subset) (JNDs)") +
  theme(axis.title = element_text(size = 10)) 
png("output//rock-v-soil//LR-v-LS-lum-gMeanDist-scatter.png", 
    width = 800, height = 600, res = 120); p; dev.off()

p <- ggplot(allGmeanDist, aes(x = dS_a_c, y = dS_a_b)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + coord_fixed(xlim = c(0,2.2), ylim = c(0,2.2)) + 
  labs(y = "L-R geometric mean chromatic distance (JNDs)",
       x = "L-S geometric mean chromatic distance (JNDs)") +
  theme(axis.title = element_text(size = 10)) 
png("output//rock-v-soil//LR-v-LS-col-gMeanDist-scatter.png", 
    width = 800, height = 600, res = 120); p; dev.off()

#* Varying nearest % of dists ####
# The nearest 10 %, 50% and all
#** Plotting data

data <- list(jnds_b, jnds_b, jnds_c, jnds_c, jnds_b, jnds_b, jnds_c, jnds_c)
locOn <- rep(T, 8)
acrIm <- rep(F,8) 
vf <- c(rep("dL",4), rep("dS", 4))
qt <- rep(c(10,2), 4)
distSubsets <- Map(nearSubset,data,locOn,acrIm, vf, qt)
names(distSubsets) <- c("b_10_dL","b_50_dL", "c_10_dL", "c_50_dL",
                        "b_10_dS", "b_50_dS", "c_10_dS", "c_50_dS")

rm(locOn,acrIm, vf, qt, data)

distSubsets <- lapply(distSubsets, distToOrigDF)
df2 <- df[,c(colnames(distSubsets[[1]]))]
df2 <- ungroup(df2)
distSubsets[["b_all_dL"]] <- df2[df2$substrate != "c",] # add in the 'all' dists
distSubsets[["b_all_dS"]] <- df2[df2$substrate != "c",]
distSubsets[["c_all_dL"]] <- df2[df2$substrate != "b",]
distSubsets[["c_all_dS"]] <- df2[df2$substrate != "b",]

rm(df2)

bindSubsets <- rbindlist(distSubsets, idcol = T) 
bindSubsets <- bindSubsets %>% 
  group_by(.id, abbrevs, substrate, mspec) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), mean) # each image contributes equally to gmean dist

distSubsets <- split(bindSubsets, bindSubsets$.id)
gMeanSubsets <- Map(f=gMeanDists,distSubsets,substrates = c(rep("b",6), rep("c", 6)))
#gMeanSubsets <- lapply(distSubsets, gMeanDists, substrates = "b") 
bindGmean <- rbindlist(gMeanSubsets, idcol = T)
bindGmean <- bindGmean %>% filter(comparison == "local") 
bindGmean <- separate(bindGmean,.id, 
                      into = c("substrate","pcentNrst","visfeature"))
bindGmean <- bindGmean %>%
  pivot_longer(cols = c("dS", "dL"),
               values_to = "gmDist", names_to = "visInfo") %>%
  filter(visfeature == visInfo) # drop rows where visInf and visf !=

gMeanSubsets <- gMeanSubsets %>% select(-c("patch1", "patch2")) %>% 
  pivot_wider(values_from = "gmDist", names_from = "substrate")

ggplot(gMeanSubsets %>% filter(visInfo == "dL"), aes(x =b, y = c))+
  geom_point()+
  #coord_equal()+
  theme(aspect.ratio = 1)+
  xlim(c(0,12))+ ylim(c(0,12))+
  facet_grid(vars(pcentNrst), scales = 'free')


#* Debugging ####

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


# bootstraps
lspop_spp <- split(df,df$pop_spp)
lsSubVec <- lapply(lspop_spp,`[`, , c('substrate'))

lspop_spp <- lapply(lspop_spp, function(x) {
  rownames(x) <- make.unique(x$substrate) ; x
})
lspop_spp <- lapply(lspop_spp, function(x) {
  x <- x[,which(names(x) %in% c("swMean", "mwMean", "lwMean", "lumMean"))] ; x
})
lspop_spp <- lapply(lspop_spp, function(x) {
  colnames(x)[which(names(x) %in% c("swMean", "mwMean", "lwMean", "lumMean"))] <- 
    c("s", "m", "l", "lum") ; x
})

bootRes <- vector(mode = "list", length = length(lspop_spp))
bootRes <- lapply(seq_along(bootRes), function(x) bootcoldist(lspop_spp[[x]],
                                                              by = lsSubVec[[x]], n = c(1,16,32), weber = .1, weber.achro = .1, 
                                                              weber.ref = "longest", qcatch = "Qi", achromatic = TRUE))
names(bootRes) <- names(lspop_spp)
bootResTransf <- lapply(bootRes, function(x) data.frame(x))
bootResTransf <- lapply(bootResTransf, function(df) transform(df, comparison = rownames(df)))
bootResTransf <- rbindlist(bootResTransf, idcol = "pop_spp")


data("sicalis")
sicmod <- vismodel(sicalis, visual = "avg.uv", relative = FALSE)
regions <- substr(rownames(sicmod), 6, 6)
sicdist <- bootcoldist(sicmod,
                       by = regions,
                       n = c(1, 2, 2, 4),
                       weber = 0.05
)

dataf <- gMeanDists(df, substrates = c("b","c"), combine_bg = F)




