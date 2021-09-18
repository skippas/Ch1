###### Local adaptation ##### 

# Within pop ie local match VERSUS
## Across pop
## Across pop within geology (divergent selection not strong, so should not differ significantly)
## Across pop across geology (divergent selection strong, so should differ significantly)
## Across species within geology (convergence or conservatism, so should not differ signficantly)
## Across species across geology (divergence / ecological speciation)

library(tidyverse)
library(data.table)
source("loading-cleaning.R")

# A single Lithops 
set.seed(42)
df <- rbind(df[df$substrate != "a",], df %>% filter(substrate == "a") %>%
               group_by(abbrevs, mspec) %>% sample_n(size = 1))
### DISTANCES ######
# nearest 10% each substrate OR nearest n ROIs across all images for a single plant 
# in each image
source("functions//coldist_effic.R")
source("functions//nearSubset.R")
jnds_c <- colDistEffic(bg_subs = "c")
jnds_b <- colDistEffic(bg_subs = "b")
nrstSub_dL_c <- nearSubset(jnds_c, localOnly = F, acrossImage = T, 
                        visfeature = "dL", qtiles =  10) 
nrstSub_dL_b <- nearSubset(jnds_b, localOnly = F, acrossImage = T, 
                        visfeature = "dL", qtiles =  10) 
nrstSub_dS_c <- nearSubset(jnds_c, localOnly = F, acrossImage = T, 
                        visfeature = "dS", qtiles =  10) 
nrstSub_dS_b <- nearSubset(jnds_b, localOnly = F, acrossImage = T, 
                        visfeature = "dS", qtiles =  10) 

transfDF <- function(df,dropInfo, keepInfo, subType){
  df %>% select(-c({{dropInfo}})) %>%
    mutate(subs = subType,
           visInfo = keepInfo) %>%
    rename("distance" = {{keepInfo}})
}
nrstSub_dL_b <- transfDF(nrstSub_dL_b, dropInfo = "dS", keepInfo = "dL", subType = "b")
nrstSub_dL_c <- transfDF(nrstSub_dL_c, dropInfo = "dS", keepInfo = "dL", subType = "c")
nrstSub_dS_b <- transfDF(nrstSub_dS_b,dropInfo = "dL", keepInfo = "dS", subType = "b")
nrstSub_dS_c <- transfDF(nrstSub_dS_c,dropInfo = "dL", keepInfo = "dS", subType  = "c")
nrstSub_bc <- rbind(nrstSub_dS_b, nrstSub_dS_c, nrstSub_dL_b, nrstSub_dL_c) 

rm(nrstSub_dL_b, nrstSub_dL_c, nrstSub_dS_b, nrstSub_dS_c, jnds_b, jnds_c)

### pairwise image nearest dists ####
perLithMeans <- nrstSub_bc %>% 
  group_by(visInfo,subs,abbrevs.x, abbrevs.y, mspec.x, mspec.y) %>% 
  summarise(pwImageDist = mean(distance)) 

### SET LITHOLOGY GROUPS #######
library(data.table) 
perLithMeans <- merge.data.table(
  setDT(perLithMeans),
  setDT(unique(
    df[, c('abbrevs', 'collapseSpp', 'allan_geol')])),
  by.x = 'abbrevs.x', by.y = 'abbrevs', all.x = T) 

perLithMeans <- perLithMeans %>%  
  rename(collapseSpp.x = 'collapseSpp',
         lithology.x = 'allan_geol')

perLithMeans <- merge.data.table(
  setDT(perLithMeans), setDT(unique(
    df[, c('abbrevs', 'collapseSpp', 'allan_geol')])),
  by.x = 'abbrevs.y', by.y = 'abbrevs', all.x = T) 

perLithMeans <- perLithMeans %>%  
  rename(collapseSpp.y = 'collapseSpp', lithology.y = 'allan_geol')

perLithMeans$comparison[with(perLithMeans, collapseSpp.x == collapseSpp.y &
                          lithology.x == lithology.y)] <- 'inSpp.inGeol'
perLithMeans$comparison[with(perLithMeans, collapseSpp.x == collapseSpp.y &
                               lithology.x != lithology.y)] <- 'inSpp.crossGeol'
perLithMeans$comparison[with(perLithMeans, collapseSpp.x != collapseSpp.y &
                               lithology.x == lithology.y)] <- 'crossSpp.inGeol'
perLithMeans$comparison[with(perLithMeans, collapseSpp.x != collapseSpp.y &
                               lithology.x != lithology.y)] <- 'crossSpp.crossGeol'
perLithMeans$comparison[with(perLithMeans, abbrevs.x == abbrevs.y)] <- 'inSite.inGeol'

perLithMeans$comparison <- factor(
  perLithMeans$comparison, levels = c(
    'inSite.inGeol', 'inSpp.inGeol',
    'inSpp.crossGeol', 'crossSpp.inGeol', 'crossSpp.crossGeol'))

# perLithMeans <- perLithMeans %>%
#   mutate(collapseSpp.x = df$collapseSpp[match(abbrevs.x,df$abbrevs)],
#          collapseSpp.y = df$collapseSpp[match(abbrevs.y,df$abbrevs)],
#          geology.x = df$allan_geol[match(abbrevs.x,df$abbrevs)],
#          geology.y = df$allan_geol[match(abbrevs.y,df$abbrevs)],
#          comparisonSpp = case_when(abbrevs.x == abbrevs.y ~ "inSite",
#                                    collapseSpp.x == collapseSpp.y ~ "inSpp",
#                                    collapseSpp.x != collapseSpp.y ~ "crossSpp"),
#          comparisonGeol = case_when(geology.x == geology.y ~ "inGeol",
#                                     geology.x != geology.y ~ "crossGeol"),
#          comparison = interaction(comparisonSpp, comparisonGeol, lex.order = T, drop = T),
#          comparison = fct_relevel(comparison,
#                                   "inSite.inGeol", "inSpp.inGeol", "inSpp.crossGeol")) %>%
#   group_by(collapseSpp.x) %>%
#   mutate(popReps = n_distinct(abbrevs.x)) %>% 
#   with_groups(NULL, mutate, collapseSpp.x = fct_reorder(collapseSpp.x, popReps, min, .desc = T)) %>%
#   select(-popReps)

### Foreign distances (paired and (more) independent) ####
perLithMeansHA <- perLithMeans %>% 
  group_by(collapseSpp.x, abbrevs.x, subs, visInfo, mspec.x, comparison) %>%
  summarise(dPerLith = mean(pwImageDist))

perLithMeansLF <- perLithMeans %>% 
  group_by(collapseSpp.y, abbrevs.y, subs, visInfo, mspec.y, comparison) %>%
  summarise(dPerRock = mean(pwImageDist))

### STATS #######

#* ANOVAs and custom contrasts ####

library(emmeans)
?"contrast-methods"
levels(perLithMeansHA$comparison) # first level automatically used as reference

anovaResLum <- perLithMeansHA %>%
  ungroup() %>% 
  filter(visInfo == "dL", subs == "b") %>%
  nest_by(abbrevs.x) %>% 
  mutate(model = list(lm(dPerLith ~ comparison, data = data)),
         contrasts = list(emmeans(model,specs =  trt.vs.ctrl ~ comparison))) %>%
  select(-c(model, data))
anovaResLum$contrasts2 <- sapply(anovaResLum$contrasts, "[[", 2)
anovaResLum$contrasts2 <- map(anovaResLum$contrasts2, as.data.frame)
anovaResLum <- anovaResLum %>% unnest(contrasts2) %>% select(-c(contrasts))
anovaResLum$contrast <- gsub(" - inSite.inGeol","", anovaResLum$contrast)

anovaResCol <- perLithMeansHA %>% ungroup() %>% 
  filter(visInfo == "dS", subs == "b") %>%
  nest_by(abbrevs.x) %>% 
  mutate(model = list(lm(dPerLith ~ comparison, data = data)),
         contrasts = list(emmeans(model,specs =  trt.vs.ctrl ~ comparison))) %>%
  select(-c(model, data))
anovaResCol$contrasts2 <-  sapply(anovaResCol$contrasts, "[[", 2)
anovaResCol$contrasts2 <-  map(anovaResCol$contrasts2, as.data.frame)
anovaResCol <- anovaResCol %>% unnest(contrasts2) %>% select(-c(contrasts))
anovaResCol$contrast <- gsub(" - inSite.inGeol","", anovaResCol$contrast)

#** Anno df ####
anovaResCol$visInfo <- "dS"
anovaResLum$visInfo <- "dL"
labelDF <- rbind(anovaResCol, anovaResLum) %>%
  rename(comparison = contrast) 
labelDF <- merge(labelDF, unique(df[c("abbrevs", "collapseSpp")]),
                 by.x = "abbrevs.x", by.y = "abbrevs", all.x = T)
labelDF$label <- "*"

labelDF <- labelDF %>% filter(p.value < 0.05) %>%
  rename(collapseSpp.x = collapseSpp) %>%
  mutate(xPos = 
           case_when(visInfo == "dS" ~ 1.5,
                     TRUE ~ 7.5)) %>%
  mutate(label = "*")


#** Tables ####
anovaResLum %>% mutate_if(is.numeric, round, 2) %>% 
  mutate(p.value = format(p.value, nsmall = 2)) %>% 
  select(c(abbrevs.x, contrast, p.value)) %>%
  pivot_wider(values_from = p.value, names_from = contrast) %>% write.csv(quote = F)
anovaResCol %>% mutate_if(is.numeric, round, 2) %>% 
  mutate(p.value = format(p.value, nsmall = 2)) %>% 
  select(c(abbrevs.x, contrast, p.value)) %>%
  pivot_wider(values_from = p.value, names_from = contrast) %>% write.csv(quote = F)

rm(anovaResCol, anovaResLum)
#* One sided T-tests, Paired / unpaired ####

perLithMeansHA <- ungroup(perLithMeansHA) %>% 
  mutate(comparison = as.character(comparison)) # fn below doesn't like the factor

ttPaired <- function(x){
  contrasts <- grep(c('Spp.cross|Spp.in'), names(x), value = T) 
  sapply(contrasts, function(y){
    t.test(x[,'inSite.inGeol'],
           x[,y], 
           alternative = 'less')$p.value})
}

pairTest <- perLithMeansHA %>% 
  pivot_wider(names_from = comparison, values_from = dPerLith) %>%
  split(.,list(.$abbrevs.x,
               .$visInfo,
               .$subs), sep = '__')

ttPairedRes <- pairTest %>%
  lapply(., function(x) 
    x[ , apply(x, 2, function(y) !any(is.na(y)))]) %>% # drop NA columns
  lapply(., ttPaired) %>%
  lapply(., as.data.frame) %>%
  lapply(., setnames, 'pVal') %>%
  lapply(., function(x){ x$pAdj <- 
    p.adjust(x$pVal, method = 'bonferroni'); x}) %>%
  lapply(., setDT, keep.rownames = 'comparison') %>%
  rbindlist(., idcol = 'abbrevs.x') %>%
  separate(., abbrevs.x, sep = '__', into = c('abbrevs.x', 'visInfo', 'subs'))

#** summary stats ####

sumStatsTt <- merge(ttPairedRes, unique(df[, c('abbrevs', 'collapseSpp')]),
                    all.x = T, by.x = 'abbrevs.x', by.y = 'abbrevs') 

# number of comparisons between different lithology categories within spp
sumStatsTt %>% 
  filter(subs == 'b', visInfo == 'dL',
         collapseSpp != 'singlePops',
         comparison != 'crossSpp.inGeol' & 
           comparison != 'crossSpp.crossGeol') %>% 
  group_by(comparison) %>%
  summarise(n = n())

# number of sig diffs to foreign lithology categories within spp
sumStatsTt %>% 
  filter(subs == 'b', visInfo == 'dL',
         collapseSpp != 'singlePops',
         comparison != 'crossSpp.inGeol' &
           comparison != 'crossSpp.crossGeol',
         pAdj < .05) %>% 
  group_by(comparison) %>%
  summarise(n = n())

sumStatsTt %>% 
  filter(collapseSpp != 'singlePops',
         subs == 'b', visInfo == 'dS',
         comparison != 'crossSpp.inGeol' & comparison != 'crossSpp.crossGeol',
         pAdj < .05) %>% 
  group_by(comparison) %>%
  summarise(n = n())

# number of across spp comparisons between lith categories
sumStatsTt %>% 
  filter(subs == 'b', visInfo == 'dL',
         comparison != 'inSpp.inGeol' & 
           comparison != 'inSpp.crossGeol') %>% 
  group_by(comparison) %>%
  summarise(n = n())

# number of sig diffs to foreign lithology categories ACROSS spp
sumStatsTt %>% 
  filter(subs == 'b', visInfo == 'dL',
         comparison != 'inSpp.inGeol' & comparison != 'inSpp.crossGeol',
         pAdj < .05) %>% 
  group_by(comparison) %>%
  summarise(n = n())

sumStatsTt %>% 
  filter(subs == 'b', visInfo == 'dS',
         comparison != 'inSpp.inGeol' & comparison != 'inSpp.crossGeol',
         pAdj < .05) %>% 
  group_by(comparison) %>%
  summarise(n = n())

# number of sig diffs to within species foreign lithologies
# WHEN both lithology categories present
ttPairedRes %>% 
  group_by(abbrevs.x) %>% 
  filter(n_distinct(comparison) > 3) %>% # all categories must be present
  filter(subs == 'b', visInfo == 'dL',
         comparison != 'crossSpp.inGeol' & comparison != 'crossSpp.crossGeol',
         pAdj < .05) %>% 
  group_by(comparison) %>%
  summarise(n = n())
  
sumStats <- perLithMeansHA %>% 
  filter(subs == 'b') %>%
  group_by(abbrevs.x, comparison, visInfo) %>%
  summarise(meanCategDist = mean(dPerLith))

sumStats %>% # which category most frequently has the lowest distance
  filter(visInfo == 'dL') %>%
  group_by(abbrevs.x) %>%
  slice_min(meanCategDist) %>%
  group_by(comparison) %>%
  summarise(n =n())

sumStats %>% # how frequently is inGeol < crossGeol, within Spp.
  group_by(abbrevs.x) %>%
  filter(n_distinct(comparison) > 4) %>%
  filter(visInfo == 'dL',
         comparison == 'inSpp.inGeol' | comparison == 'inSpp.crossGeol') %>%
  slice_min(meanCategDist) %>%
  group_by(comparison) %>%
  summarise(n =n())

sumStats %>% # how frequently is inGeol < crossGeol, across Spp.
  group_by(abbrevs.x) %>%
  filter(visInfo == 'dL',
         comparison == 'crossSpp.inGeol' | 
           comparison == 'crossSpp.crossGeol') %>%
  slice_min(meanCategDist) %>%
  group_by(comparison) %>%
  summarise(n =n())

rm(sumStats)
#** Anno df ####

# doesn't the below cause multiple label issues (label repeat for each row?)
perLithMeansHA <- merge(perLithMeansHA, ttPairedRes, all.x = T) %>% 
  mutate(xPos = 
           case_when(visInfo == "dS" ~ 1.8,
                     TRUE ~ 10),
         label = case_when(pAdj < .05 ~ '*',
                           TRUE ~ ''))
rm(ttPaired, pairTest)

### PLOTS ##########


# Plot setup 
theme_set(theme_bw())
theme_update(panel.spacing = unit(0, "lines"))

library(RColorBrewer)
library(scales)
show_col(palette.colors(palette = 'okabe-ito'))

perLithMeansHA$comparison <- factor(
  perLithMeansHA$comparison, levels = c(
    'inSite.inGeol', 'inSpp.inGeol',
    'inSpp.crossGeol', 'crossSpp.inGeol', 'crossSpp.crossGeol'))

myColors <- palette.colors(palette = 'okabe-ito')[c(1,7,3,7,3)]
names(myColors) <- levels(perLithMeansHA$comparison)
myShapes <- c(16,16,16,15,15)
names(myShapes) <- levels(perLithMeansHA$comparison)

winSppPlots <- c('inSite.inGeol', 'inSpp.inGeol', 'inSpp.crossGeol')
crossSppPlots <- c('inSite.inGeol','crossSpp.inGeol','crossSpp.crossGeol')
visInfoRename <- c(dL = 'Luminance', dS = 'Chroma')
perLithMeansHA[perLithMeansHA$subs == 'b','subRename'] <- 'rock'
perLithMeansHA[perLithMeansHA$subs == 'c', 'subRename'] <- 'soil'
perLithMeansHA <- perLithMeansHA[
  with(perLithMeansHA, !(collapseSpp.x == 'singlePops' &
                           (comparison == 'inSpp.inGeol' |
                         comparison == 'inSpp.crossGeol'))),]
perLithMeansHA <- merge(
  perLithMeansHA, unique(df[,c('abbrevs','geol_abbrevs')]),
      by.x = 'abbrevs.x', by.y = 'abbrevs', all.x = T) %>%
  separate(., 'abbrevs.x', into = 'abbrevs.x', sep = '_') %>%
  mutate(abbrevs.x = paste(.$abbrevs.x, .$geol_abbrevs, sep = '_'))

pd1 = position_dodge(0.5)

plotFn <-function(df, withinOrAcr, RockSoil){
  p <- ggplot(df %>%
                filter(comparison %in% withinOrAcr, subRename == RockSoil),
              aes(x = dPerLith, y = abbrevs.x, color=comparison,
                  shape = comparison))+
    stat_summary(fun.data= mean_cl_boot, geom="errorbar", 
                 position=pd1, width = 0.5)+
    stat_summary(fun = mean, geom="point", size=1, position=pd1)+
    geom_text(aes(x = xPos, label = label),
              position = pd1, show.legend = F)+
    facet_grid(collapseSpp.x~visInfo,
               scales = "free", space = "free_y",
               labeller = labeller(visInfo = visInfoRename))+
    xlab(paste(RockSoil, 'distance (JNDs)'))+
    ylab("Lithops populations")+
    scale_shape_manual(name = paste(RockSoil, 'groups'),
                       values = myShapes,
                       limits = force)+
    scale_colour_manual(name = paste(RockSoil,'groups'),
                        values = myColors,
                        limits = force)
}
p <- plotFn(perLithMeansHA, winSppPlots, 'rock')
png('output//local-adaptation//withinSpp_loc-v-forLithology_rock.png',
    type = 'cairo', width = 20, height = 30, units = 'cm', res = 300,
    pointsize = 6); p ; dev.off()
p <- plotFn(perLithMeansHA, crossSppPlots, 'rock')
png('output//local-adaptation//crossSpp_loc-v-forLithology_rock.png',
    type = 'cairo', width = 20, height = 30, units = 'cm', res = 300,
    pointsize = 6); p ; dev.off()

p <- plotFn(perLithMeansHA, winSppPlots, 'soil')
png('output//local-adaptation//withinSpp_loc-v-forLithology_soil.png',
    type = 'cairo', width = 20, height = 30, units = 'cm', res = 300,
    pointsize = 6); p ; dev.off()
p <- plotFn(perLithMeansHA, crossSppPlots, 'soil')
png('output//local-adaptation//crossSpp_loc-v-forLithology_soil.png',
    type = 'cairo', width = 20, height = 30, units = 'cm', res = 300,
    pointsize = 6); p ; dev.off()
rm(p, pd1)


### LF version ####
# STATS
#* One sided T-tests, Paired / unpaired 

perLithMeansLF <- ungroup(perLithMeansLF) %>% 
  mutate(comparison = as.character(comparison)) # fn below doesn't like the factor

ttPaired <- function(x){
  contrasts <- grep(c('Spp.cross|Spp.in'), names(x), value = T) 
  sapply(contrasts, function(y){
    t.test(x[,'inSite.inGeol'],
           x[,y], 
           alternative = 'less')$p.value})
}

pairTest <- perLithMeansLF %>% 
  pivot_wider(names_from = comparison, values_from = dPerRock) %>%
  split(.,list(.$abbrevs.y,
               .$visInfo,
               .$subs), sep = '__')

ttPairedRes <- pairTest %>%
  lapply(., function(x) 
    x[ , apply(x, 2, function(y) !any(is.na(y)))]) %>% # drop NA columns
  lapply(., ttPaired) %>%
  lapply(., as.data.frame) %>%
  lapply(., setnames, 'pVal') %>%
  lapply(., function(x){ x$pAdj <- 
    p.adjust(x$pVal, method = 'bonferroni'); x}) %>%
  lapply(., setDT, keep.rownames = 'comparison') %>%
  rbindlist(., idcol = 'abbrevs.y') %>%
  separate(., abbrevs.y, sep = '__', into = c('abbrevs.y', 'visInfo', 'subs'))
# Anno df 

# doesn't the below cause multiple label issues (label repeat for each row?)
perLithMeansLF <- merge(perLithMeansLF, ttPairedRes, all.x = T) %>% 
  mutate(xPos = 
           case_when(visInfo == "dS" ~ 1.8,
                     TRUE ~ 10),
         label = case_when(pAdj < .05 ~ '*',
                           TRUE ~ ''))
rm(ttPaired, pairTest)

# PLOTS 

# Plot setup 
theme_set(theme_bw())
theme_update(panel.spacing = unit(0, "lines"))

library(RColorBrewer)
library(scales)
show_col(palette.colors(palette = 'okabe-ito'))

perLithMeansLF$comparison <- factor(
  perLithMeansLF$comparison, levels = c(
    'inSite.inGeol', 'inSpp.inGeol',
    'inSpp.crossGeol', 'crossSpp.inGeol', 'crossSpp.crossGeol'))

myColors <- palette.colors(palette = 'okabe-ito')[c(1,7,3,7,3)]
names(myColors) <- levels(perLithMeansLF$comparison)
myShapes <- c(16,16,16,15,15)
names(myShapes) <- levels(perLithMeansLF$comparison)

winSppPlots <- c('inSite.inGeol', 'inSpp.inGeol', 'inSpp.crossGeol')
crossSppPlots <- c('inSite.inGeol','crossSpp.inGeol','crossSpp.crossGeol')
visInfoRename <- c(dL = 'Luminance', dS = 'Chroma')
perLithMeansLF[perLithMeansLF$subs == 'b','subRename'] <- 'rock'
perLithMeansLF[perLithMeansLF$subs == 'c', 'subRename'] <- 'soil'
perLithMeansLF <- perLithMeansLF[
  with(perLithMeansLF, !(collapseSpp.y == 'singlePops' &
                           (comparison == 'inSpp.inGeol' |
                              comparison == 'inSpp.crossGeol'))),]
perLithMeansLF <- merge(
  perLithMeansLF, unique(df[,c('abbrevs','geol_abbrevs')]),
  by.x = 'abbrevs.y', by.y = 'abbrevs', all.x = T) %>%
  separate(., 'abbrevs.y', into = 'abbrevs.y', sep = '_') %>%
  mutate(abbrevs.y = paste(.$abbrevs.y, .$geol_abbrevs, sep = '_'))

pd1 = position_dodge(0.5)

plotFn <-function(df, withinOrAcr, RockSoil){
  p <- ggplot(df %>%
                filter(comparison %in% withinOrAcr, subRename == RockSoil),
              aes(x = dPerRock, y = abbrevs.y, color=comparison,
                  shape = comparison))+
    stat_summary(fun.data= mean_cl_boot, geom="errorbar", 
                 position=pd1, width = 0.5)+
    stat_summary(fun = mean, geom="point", size=1, position=pd1)+
    geom_text(aes(x = xPos, label = label),
              position = pd1, show.legend = F)+
    facet_grid(collapseSpp.y~visInfo,
               scales = "free", space = "free_y",
               labeller = labeller(visInfo = visInfoRename))+
    xlab(paste('Distance to Lithops groups (JNDs)'))+
    ylab(paste(RockSoil,"populations"))+
    scale_shape_manual(name = paste('Lithops groups'),
                       values = myShapes,
                       limits = force)+
    scale_colour_manual(name = paste('Lithops groups'),
                        values = myColors,
                        limits = force)
}
p <- plotFn(perLithMeansLF, winSppPlots, 'rock')
png('output//local-adaptation//LF_withinSpp_loc-v-forLithology_rock.png',
    type = 'cairo', width = 20, height = 30, units = 'cm', res = 300,
    pointsize = 6); p ; dev.off()
p <- plotFn(perLithMeansLF, crossSppPlots, 'rock')
png('output//local-adaptation//LF_crossSpp_loc-v-forLithology_rock.png',
    type = 'cairo', width = 20, height = 30, units = 'cm', res = 300,
    pointsize = 6); p ; dev.off()

p <- plotFn(perLithMeansLF, winSppPlots, 'soil')
png('output//local-adaptation//LF_withinSpp_loc-v-forLithology_soil.png',
    type = 'cairo', width = 20, height = 30, units = 'cm', res = 300,
    pointsize = 6); p ; dev.off()
p <- plotFn(perLithMeansLF, crossSppPlots, 'soil')
png('output//local-adaptation//LF_crossSpp_loc-v-forLithology_soil.png',
    type = 'cairo', width = 20, height = 30, units = 'cm', res = 300,
    pointsize = 6); p ; dev.off()
rm(p, pd1)


### Overall Ttests background categories  ################

#* Setup plotting dataframes ####
transfDF <- function(visInf, substr){
  
  overallContrasts <- perLithMeansHA %>% 
    filter(subRename == substr, visInfo == visInf) %>% 
    group_by(abbrevs.x, comparison) %>%
    summarise(meanDist = mean(dPerLith)) %>%
    mutate(comparison = as.factor(comparison)) %>%
    mutate(comparison = fct_relevel(comparison, "inSite.inGeol", "inSpp.inGeol", 
                                    "inSpp.crossGeol", "crossSpp.inGeol")) 
}
lumContr <- transfDF("dL","rock") # How do NA values affect the tests?
colContr <- transfDF("dS", "rock")

#* Run paired t-tests ####
custContrasts <- function(contrastDF){
  
  ttestCompat <- contrastDF %>% 
    pivot_wider(names_from = comparison, values_from = meanDist)
  contrasts <- c("inSpp.inGeol", "inSpp.crossGeol",
                 "crossSpp.crossGeol","crossSpp.inGeol")   
  tResults <- lapply(contrasts, function(x){
    t.test(ttestCompat$inSite.inGeol,
           ttestCompat[[x]], paired = T)
  })
  names(tResults) <- contrasts
  return(tResults)
}
lumRes <- custContrasts(lumContr)
colRes <- custContrasts(colContr)


#* annotation df ####
statLabeller <- function(tResults, yPos){
  
  statLabels <- data.frame(.y. = "meanDist", group1 = "inSite.inGeol",
                           group2 = c(names(tResults)),
                           p = unname(unlist(lapply(tResults,
                                                    `[[`, "p.value"))),
                           statistic = unname(unlist(lapply(tResults,
                                                    `[[`, "statistic"))),
                           df = unname(unlist(lapply(tResults,
                                                     `[[`, "parameter"))))
  statLabels$p.adj <- p.adjust(statLabels$p, n = length(tResults),
                               method = "bonferroni")
  statLabels$p.format <- format.pval(statLabels$p.adj,digits = 2)
  statLabels$statistic <- round(statLabels$statistic, digits = 2)
  statLabels$label <- with(statLabels, paste('p =',p.format,'  ',
                            't =',statistic,'  ', 'df =',df ))
  statLabels <- statLabels %>%
    mutate(y.position = case_when(group2 == "inSpp.inGeol" ~ yPos[1],
                                  group2 == "inSpp.crossGeol" ~ yPos[2],
                                  group2 == "crossSpp.inGeol" ~ yPos[3],
                                  group2 == "crossSpp.crossGeol" ~ yPos[4]))
  return(statLabels)
}
lumLabels <- statLabeller(lumRes, yPos =  c(5:8))
colLabels <- statLabeller(colRes, yPos =  seq(1,1.5, length = 4))
meanLumLabels <- lumContr %>% group_by(comparison) %>%
  summarise(meanDist = round(mean(meanDist), digits = 2))
meanChrLabels <- colContr %>% group_by(comparison) %>%
  summarise(meanDist = round(mean(meanDist), digits = 2))

#* overall plots ####
library(ggpubr)

nCateg <- lumContr %>% ungroup() %>% count(comparison)
bgCategLabels <- c(inSite.inGeol = 'local', 
                   inSpp.inGeol = 'same species & lithology',
                   inSpp.crossGeol = 'same species, \ndifferent lithology',
                   crossSpp.inGeol = 'different species, \nsame lithology',
                   crossSpp.crossGeol = 'different species & lithology')

pLum <- ggboxplot(lumContr, x = "comparison", y = "meanDist", width = 0.3,
               outlier.shape = 1)+
  stat_pvalue_manual(lumLabels, label = 'label', size = 2.5)+
  geom_text(data = meanLumLabels,
            aes(y = meanDist, label = meanDist), size = 2)+
  scale_y_continuous(limits = c(0,9))+
  ylab("Mean luminance distances (JNDs) \nto background categories")+
  xlab("")+ 
  theme_minimal()

pChr <- ggboxplot(colContr, x = "comparison", y = "meanDist", width = 0.3,
               outlier.shape = 1)+
  stat_pvalue_manual(colLabels, label = 'label', size = 2.5)+ 
  geom_text(data = meanChrLabels,
            aes(y = meanDist-0.008, label = meanDist), size = 2)+
  ylab("Mean chromatic distances (JNDs) \nto background categories")+
  xlab("Species x lithology \nbackground categories")+
  theme_minimal()

png("output//local-adaptation//Overall-bgCategories-Ttests-lumChr.png", type = 'cairo',
    width = 215, height = 140, units = "mm", res = 300, pointsize = 6)
pLum + pChr + plot_layout(nrow = 2)
dev.off()

#* across all populations (old figures) ####
pdata <- perLithMeans %>% filter(visInfo == "dL") %>%
  group_by(abbrevs.y,abbrevs.x) %>% 
  summarise(meanDPerLith = mean(dPerLith))
ggplot() + 
  geom_jitter(data = pdata %>% filter(abbrevs.x != abbrevs.y),
              aes(y = abbrevs.y, x = meanDPerLith)) +
  geom_jitter(data = pdata %>% filter(abbrevs.x == abbrevs.y), 
              aes(y = abbrevs.y, x = meanDPerLith, color = "red"))
rm(pdata)

png(file = 'output//lum-near-dist-medians.png', width = 5, height = 5, units = "cm")
ggplot(perLithMeans, aes(x = abbrevs.y,y=median)) +
  geom_jitter(data = medsLum %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = medsLum %>% filter(group == "local"),width=0.1,alpha=0.8,color = "red", size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=0.95,vjust=0.9, size = 8, face = "bold")) + 
  #theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")
dev.off()



### Industrial plots (LF, HA, RS, LUMCOL, etc. different combs ) ####

allRS <- perLithMeans
winRS <- perLithMeans %>% 
  filter(comparison != "crossSpp.crossGeol" & comparison != 'crossSpp.inGeol',
         collapseSpp.x != "singlePops",
         collapseSpp.y != "singlePops") # there are no within spp comps?
acrossRS <- perLithMeans %>%
  filter(comparison != "inSpp.inGeol" & comparison != 'inSpp.crossGeol') 

allRock <- allRS %>% filter(subRename == "rock")
winRock <- winRS %>% filter(subRename == "rock")
acrossRock <- acrossRS %>% filter(subRename =="rock")

criterionTransf <- function(df){ # (Im not sure this is correct HA LF values) store both HA and LF criteria in single df
  critHA <- df
  critLF <- df
  critHA$criterion <- "HvA"
  critLF$criterion <- "LvF"
  critLF <- critLF %>% rename(abbrevs = abbrevs.y,
                              abbrevs.y = abbrevs.x,
                              collapseSpp = collapseSpp.y,
                              collapseSpp.y = collapseSpp.x) %>%
    rename(abbrevs.x = abbrevs,
           collapseSpp.x = collapseSpp)
  critComp <- rbind(critHA, critLF)
  return(critComp)
}

HALFwin <- criterionTransf(winRock)
HALFacross <- criterionTransf(acrossRock)


#* Run industrial plots ####

# General plot function 
# HA criterion - yaxis = abbrevs.x and facets = all distances to a lithops population
pd1 = position_dodge(0.5)

plotF <- function(df,visInf, myAxisTitle,hFacet){
  
  p <- ggplot(df %>% filter(visInfo %in% visInf),
              aes(x = dPerLith, y = abbrevs.x, color=comparison))+
    stat_summary(fun.data= mean_cl_boot, geom="errorbar", 
                 position=pd1, aes(color = comparison), width = 0.5)+
    stat_summary(fun = mean, geom="point", size=1, position=pd1)+
    facet_grid(collapseSpp.x~.data[[hFacet]], scales = "free", space = "free_y")+
    colScale+ xlab(myAxisTitle)+ ylab("Populations")
  
  # Annotate with stats if applicable
  if(length(unique(df$subRename)) < 2 && !"criterion" %in% colnames(df)){
    # remove categories from labelDF if not in the plotting df
    df <- df %>% ungroup() %>% as.data.frame()
    labelz <- labelDF[labelDF$comparison %in% unique(df[,"comparison"]),]
    labelz <- labelz[labelz$collapseSpp.x %in% unique(df[,"collapseSpp.x"]),]
    p <- p +geom_text(data = labelz,aes(x = xPos, label = label), position = pd1,
                      show.legend = F)
  }
  return(p)
}


RS_dL_all <- list(allRS, "dL", myAxisTitle = "Luminance distance", hFacet = "subRename")
RS_dS_all <- list(allRS, "dS", myAxisTitle = "chromatic distance", hFacet = "subRename")
RS_dS_win <- list(winRS, "dS", myAxisTitle = "Chromatic distance", hFacet = "subRename")
RS_dL_win <- list(winRS, "dL", myAxisTitle = "Luminance distance", hFacet = "subRename")
RS_dS_across <- list(acrossRS, "dS", myAxisTitle = "Chromatic distance", hFacet = "subRename")
RS_dL_across <- list(acrossRS, "dL", myAxisTitle = "Luminance distance", hFacet = "subRename")
lumCol_win <- list(winRock, c("dS","dL"), myAxisTitle = "Distance", hFacet = "visInfo")
lumCol_across <- list(acrossRock, c("dS", "dL"), myAxisTitle = "Distance", hFacet = "visInfo")
lumCol_all <- list(allRock, c("dS","dL"), myAxisTitle = "Distance", hFacet = "visInfo")
HALF_dL_win <- list(HALFwin, "dL", myAxisTitle = "Luminance distance", hFacet = "criterion")
HALF_dS_win <- list(HALFwin, "dS", myAxisTitle = "Chromatic distance", hFacet = "criterion")
HALF_dL_across <- list(HALFacross, "dL", myAxisTitle = "Luminance distance", hFacet = "criterion")
HALF_dS_across <- list(HALFacross, "dS", myAxisTitle = "Chromatic distance", hFacet = "criterion")

patternz <- c("^RS_", "^lumCol_", "HALF_")
lsPlots <- unlist(lapply(patternz, function(x)
  ls(pattern = x,envir=parent.env(environment()))))
lsPlots <- mget(lsPlots)

myplots <- vector('list', length(lsPlots))
names(myplots) <- names(lsPlots)
for(i in 1:length(lsPlots)){
  myplots[[i]] <- plotF(lsPlots[[i]][[1]],lsPlots[[i]][[2]],
                        lsPlots[[i]][[3]], lsPlots[[i]][[4]])
  png(paste0('output//local-adaptation//',names(myplots[i]), '.png'),
      type = 'cairo', width = 20, height = 30, units = 'cm', res = 300,
      pointsize = 6)
  print(myplots[[i]])
  dev.off()
}

patternz <- c(patternz, "^HALF", "^across", "^win","^label", "^anovaRes", "^all")
bin <- lapply(patternz, function(x)
  ls(pattern = x, envir=parent.env(environment())))
rm(list = unlist(bin), bin, patternz, myplots, lsPlots,pd1, colScale)

