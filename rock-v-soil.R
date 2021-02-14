rm(list = ls())
library(gridExtra)
library(data.table)
library(pavo)

source("functions//coldist_effic.R")
source("functions//nearest_qtile.R")

jnds <- colDistEffic(bg_subs = "b")
medsLum_b <- nearestQtile(df = jnds,visfeature = "dL", ntiles = 4, nearest_ntile = 1) 
jnds <- colDistEffic(bg_subs = "c")
medsLum_c <- nearestQtile(df = jnds,visfeature = "dL", ntiles = 4, nearest_ntile = 1) 

# brightness 

locMedsLum_b <- medsLum_b[medsLum_b$abbrevs.y == medsLum_b$abbrevs.x,]
locMedsLum_c <- medsLum_c[medsLum_c$abbrevs.y == medsLum_c$abbrevs.x,]
locMedsLum_b$median_b <- locMedsLum_b$median
locMedsLum_c$median_c <- locMedsLum_c$median
locMedsLum_bc <- merge(locMedsLum_b, by.x = "abbrevs.x", by.y = "abbrevs.x", locMedsLum_c, all.y = TRUE)

texts <- locMedsLum_bc %>% filter(median_b > 1 | median_c > 5)
png(file = 'output//rock-v-soil//brightness-match-rock-v-soil.png', width = 1024, height = 768, units = "px")
ggplot(locMedsLum_bc, aes(x = median_b, y = median_c)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  #geom_text(data = texts, aes(label = abbrevs.x,y = median_c, x = median_b)) +
  theme_bw() + coord_fixed(xlim = c(0,40), ylim = c(0,40)) + 
  labs(y = "achro dists to soil", x = "achro dists to rock") +
  theme(axis.title = element_text(size = 20)) 
dev.off()

# colour

jnds <- colDistEffic(bg_subs = "b")
medsCol_b <- nearestQtile(df = jnds,visfeature = "dS", ntiles = 4, nearest_ntile = 1) 
jnds <- colDistEffic(bg_subs = "c")
medsCol_c <- nearestQtile(df = jnds,visfeature = "dS", ntiles = 4, nearest_ntile = 1) 


locMedsCol_b <- medsCol_b[medsCol_b$abbrevs.y == medsCol_b$abbrevs.x,]
locMedsCol_c <- medsCol_c[medsCol_c$abbrevs.y == medsCol_c$abbrevs.x,]
locMedsCol_b$median_b <- locMedsCol_b$median
locMedsCol_c$median_c <- locMedsCol_c$median
locMedsCol_bc <- merge(locMedsCol_b, by.x = "abbrevs.x",
                       by.y = "abbrevs.x", locMedsCol_c, all.y = TRUE)

texts <- locMedsCol_bc %>% filter(median_b > 1 | median_c > 5)
png(file = 'output//rock-v-soil//color-match-rock-v-soil.png', width = 1024, height = 768, units = "px")
ggplot(locMedsCol_bc, aes(x = median_b, y = median_c)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  #geom_text(data = texts, aes(label = abbrevs.x,y = median_c, x = median_b)) +
  theme_bw() + coord_fixed(xlim = c(0,9), ylim = c(0,9)) + 
  labs(y = "chrom dists to soil", x = "chrom dists to rock") +
  theme(axis.title = element_text(size = 20)) 
dev.off()
