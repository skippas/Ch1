# order populations by:
# 1. distance between the local bg substrates
# 2. by the spread of substrates
# 3. by species
# 4. by geology
source("loading-cleaning.r")
library(tidyverse)
library("plotrix")

df$logLum <- exp(log(df$lumMean)) # whats the math effect of doing this!!!
sumStats <- df %>% group_by(pop_spp) %>% filter(substrate == "b") %>%
  summarise(rock_sd = sd(logLum), rock_n = n())

grpDF <- df %>% group_by(abbrevs, substrate)

# order by rock soil mean dist (var in adaptive landscape)
subLumMeans <- summarise(grpDF,logLumMean = mean(logLum))
subLumMeans <- pivot_wider(subLumMeans, names_from = substrate,
                           values_from = logLumMean)
subLumMeans$lumDist <- subLumMeans$b - subLumMeans$c
subLumMeans$abbrevs <- fct_reorder(subLumMeans$abbrevs, subLumMeans$lumDist, min)
df$abbrevs <- fct_relevel(df$abbrevs, levels(subLumMeans$abbrevs))

# order by spread of rocks
source("functions//MAD.R")
spread <- grpDF %>% filter(substrate == "b") %>% summarise(lumMAD = MAD(logLum))
spread$abbrevs <- fct_reorder(spread$abbrevs, spread$lumMAD)
df$abbrevs <- fct_relevel(df$abbrevs, levels(spread$abbrevs))
