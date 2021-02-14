source("loading-cleaning.R")
# rock-soil dep rel amounts
# brightness
rm(list = ls())
source("functions//coldist_effic.R")
source("functions//nearest_qtile.R")

jnds <- colDistEffic(bg_subs = "b")
medsLumB <- nearestQtile(df = jnds,visfeature = "dL",ntiles = 4, nearest_ntile = 1)
jnds <- colDistEffic(bg_subs = "c")
medsLumC <- nearestQtile(df = jnds,visfeature = "dL",ntiles = 4, nearest_ntile = 1)
colnames(medsLumB)[which(names(medsLumB) == "median")] <- "median_b"
colnames(medsLumC)[which(names(medsLumC) == "median")] <- "median_c"
medsLum <- merge(medsLumB,medsLumC)
medsLum$group <- "foreign"
medsLum[medsLum$abbrevs.x ==  medsLum$abbrevs.y,"group"] <- "local"
medsLum <- medsLum[medsLum$abbrevs.x ==  medsLum$abbrevs.y,]

source("rock-estimate.R")

props <- df %>% count(abbrevs, substrate)
props <- props %>% pivot_wider(names_from = substrate, values_from = n)
props$rock_cover <- (props$b / (props$b + props$c)*100)
props$soil_cover <- (props$c / (props$b + props$c)*100)
props <- arrange(props, as.factor(abbrevs)) 

propDist <- merge(props, medsLum[,c(2,3,4)], by.x = "abbrevs", by.y = "abbrevs.x")

png(file = 'output//rock-v-soil//Lum_c_dep-surface-area.png', width = 1024, height = 768, units = "px")
ggplot(propDist, aes(x = soil_cover, y = median_c)) +
  geom_point() + theme_bw() + geom_smooth(method = "lm") + 
  ylab("achro dists to soil") + xlab("% of surface area = soil")
dev.off()

png(file = 'output//rock-v-soil//Lum_b_dep-surface-area.png', width = 1024, height = 768, units = "px")
ggplot(propDist, aes(x = rock_cover, y = median_b)) +
  geom_point() + theme_bw() + geom_smooth(method = "lm") + 
  ylab("achro dists to rock") + xlab("% of surface area = rock")
dev.off()

summary(lm(median_c~soil_cover, data = propDist))
summary(lm(median_b~rock_cover, data = propDist))

library(ggrepel)
texts <- propDist[propDist$median > 1,]
ggplot(propDist,aes(x = rock_cover, y = median, label = abbrevs)) +
  geom_point() +
  geom_text_repel(data = subset(propDist, median > 1)) +
  theme_bw() + geom_smooth(method = "lm") 

  
# color

rm(list = ls())
source("functions//coldist_effic.R")
source("functions//nearest_qtile.R")

jnds <- colDistEffic(bg_subs = "b")
medsLumB <- nearestQtile(df = jnds,visfeature = "dS",ntiles = 4, nearest_ntile = 1)
jnds <- colDistEffic(bg_subs = "c")
medsLumC <- nearestQtile(df = jnds,visfeature = "dS",ntiles = 4, nearest_ntile = 1)
colnames(medsLumB)[which(names(medsLumB) == "median")] <- "median_b"
colnames(medsLumC)[which(names(medsLumC) == "median")] <- "median_c"
medsLum <- merge(medsLumB,medsLumC)
medsLum$group <- "foreign"
medsLum[medsLum$abbrevs.x ==  medsLum$abbrevs.y,"group"] <- "local"
medsLum <- medsLum[medsLum$abbrevs.x ==  medsLum$abbrevs.y,]

source("rock-estimate.R")

props <- df %>% count(abbrevs, substrate)
props <- props %>% pivot_wider(names_from = substrate, values_from = n)
props$rock_cover <- (props$b / (props$b + props$c)*100)
props$soil_cover <- (props$c / (props$b + props$c)*100)
props <- arrange(props, as.factor(abbrevs)) 

propDist <- merge(props, medsLum[,c(2,3,4)], by.x = "abbrevs", by.y = "abbrevs.x")

png(file = 'output//rock-v-soil//Col_c_dep-surface-area.png', width = 1024, height = 768, units = "px")
ggplot(propDist, aes(x = rock_cover, y = median_b)) +
  geom_point() + theme_bw() + geom_smooth(method = "lm") + 
  ylab("chrom dists to rock") + xlab("% of surface area = rock")
dev.off()

png(file = 'output//rock-v-soil//Col_b_dep-surface-area.png', width = 1024, height = 768, units = "px")
ggplot(propDist, aes(x = soil_cover, y = median_c)) +
  geom_point() + theme_bw() + geom_smooth(method = "lm") + 
  ylab("chrom dists to soil") + xlab("% of surface area = soil")
dev.off()

summary(lm(median_c~soil_cover, data = propDist))
summary(lm(median_b~rock_cover, data = propDist))

library(ggrepel)
texts <- propDist[propDist$median > 1,]
ggplot(propDist,aes(x = rock_cover, y = median, label = abbrevs)) +
  geom_point() +
  geom_text_repel(data = subset(propDist, median > 1)) +
  theme_bw() + geom_smooth(method = "lm") 



