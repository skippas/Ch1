rm(list = ls())
library(rlang)
source("functions//coldist_effic.R")
source("functions//nearest_qtile.R")

jnds <- colDistEffic(bg_subs = "b")

color_b <- vector("list",4)
lum_b <- vector("list",4)
for(i in 1:4) {
color_b[[i]] <- nearestQtile(df = jnds,visfeature = "dS",ntiles = 4, nearest_ntile = i)
lum_b[[i]] <- nearestQtile(df = jnds,visfeature = "dL",ntiles = 4, nearest_ntile = i)
}

jnds <- colDistEffic(bg_subs = "c")

color_c <- vector("list",4)
lum_c <- vector("list",4)
for(i in 1:4) {
  color_c[[i]] <- nearestQtile(df = jnds,visfeature = "dS",ntiles = 4, nearest_ntile = i)
  lum_c[[i]] <- nearestQtile(df = jnds,visfeature = "dL",ntiles = 4, nearest_ntile = i)
}

plots_col <- vector("list",4)
plots_lum <- vector("list",4)
xlabel <- c("median JND to rock", "", "", "")
ylabel <- c("median JND to soil", "", "", "")
x <- c(1/4*100, 2/4*100, 3/4*100, 4/4*100)
anno <- c("Nearest ", "x", "% of distances")
for(i in 1:length(color_c)){
  # color
  anno[2] <- x[i]
  MedsCol_c <- color_c[[i]]
  MedsCol_c <- MedsCol_c[MedsCol_c$abbrevs.y == MedsCol_c$abbrevs.x,]
  
  MedsCol_b <- color_b[[i]]
  MedsCol_b <- MedsCol_b[MedsCol_b$abbrevs.y == MedsCol_b$abbrevs.x,]
  
  MedsCol_c$median_c <- MedsCol_c$median
  MedsCol_b$median_b <- MedsCol_b$median
  
  MedsCol_bc <- merge(MedsCol_b, MedsCol_c[, c("abbrevs.x", "median_c")],by.x = "abbrevs.x",
                         by.y = "abbrevs.x", all.y = TRUE)

  plots_col[[i]] <- ggplot(MedsCol_bc, aes(x = median_b, y = median_c)) + geom_point() +
    #geom_text(data = texts, aes(label = abbrevs.x,y = median_c, x = median_b)) +
    theme_bw() + coord_fixed() + xlim(c(0,15)) + ylim(c(0,15)) + 
    labs(y = ylabel[i], x = xlabel[i]) +
    geom_hline(yintercept = 1, linetype = 'dotted', col = 'red') +
    geom_vline(xintercept = 1, linetype = 'dotted', col = 'red') +
    geom_abline(slope = 1,linetype = 'dashed') +
    annotate("text", y= 14, x =5,label= paste(anno, collapse = ""), size = 3)

  #brightness
  MedsLum_c <- lum_c[[i]]
  MedsLum_c <- MedsLum_c[MedsLum_c$abbrevs.y == MedsLum_c$abbrevs.x,]
  
  MedsLum_b <- lum_b[[i]]
  MedsLum_b <- MedsLum_b[MedsLum_b$abbrevs.y == MedsLum_b$abbrevs.x,]
  
  MedsLum_c$median_c <- MedsLum_c$median
  MedsLum_b$median_b <- MedsLum_b$median
  
  MedsLum_bc <- merge(MedsLum_b,MedsLum_c[,c("abbrevs.x","median_c")], by.x = "abbrevs.x",
                      by.y = "abbrevs.x", all.y = TRUE)
  
  plots_lum[[i]] <- ggplot(MedsLum_bc, aes(x = median_b, y = median_c)) + geom_point() +
    #geom_text(data = texts, aes(label = abbrevs.x,y = median_c, x = median_b)) +
    theme_bw() + coord_fixed() + xlim(c(0,60)) + ylim(c(0,60)) +
    labs(y = ylabel[i], x = xlabel[i]) +
    geom_hline(yintercept = 1, linetype = 'dotted', col = 'red') +
    geom_vline(xintercept = 1, linetype = 'dotted', col = 'red') +
    geom_abline(slope = 1,linetype = 'dashed') +
    annotate("text", y= 60, x =20,label= paste(anno, collapse = ""), size = 3)
}

grid.arrange(plots_lum[[1]], plots_lum[[2]],plots_lum[[3]], plots_lum[[4]], ncol=4)
grid.arrange(plots_col[[1]], plots_col[[2]],plots_col[[3]], plots_col[[4]], ncol=4)

# local v nonlocal
jnds <- colDistEffic()

medsCol <- vector("list",4)
medsLum <- vector("list",4)

for(i in 1:4) {
  medsCol[[i]] <- nearestQtile(df = jnds,visfeature = "dS",ntiles = 4, nearest_ntile = i)
  medsLum[[i]] <- nearestQtile(df = jnds,visfeature = "dL",ntiles = 4, nearest_ntile = i)
}

plots_lum <- vector("list",4)
plots_lumHist <- vector("list",4)

for(i in 1:length(medsLum)){ 
  
  lum <- medsLum[[i]]
  lum$abbrevs <- with(lum, reorder(abbrevs.y, median, function(x) max(x)))
  
  plots_lum[[i]] <- ggplot(lum, aes(x = abbrevs.x,y=median)) +
    geom_jitter(data = lum %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + 
    geom_jitter(data = lum %>% filter(group == "local"),width=0.1,alpha=0.8,color = "red", size = 2) + theme_bw() +
    theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
    labs(y = "lithops median distance to background", x = "background populations")
  plots_lumHist[[i]] <- ggplot(lum) + geom_freqpoly(aes(x = median, colour = group), stat = "density")

}
temp <- medsLum

names(medsLum) <- c(1/4*100, 2/4*100, 3/4*100, 4/4*100)
for(i in 1:length(medsLum)){
  medsLum[[i]]$qtile <- names(medsLum[i])
  medsLum[[i]] <- na.omit(medsLum[[i]])
  medsLum[[i]]$group <- "0"
  medsLum[[i]][medsLum[[i]]$abbrevs.x ==  medsLum[[i]]$abbrevs.y,"group"] <- "local"
  medsLum[[i]][medsLum[[i]]$abbrevs.x !=  medsLum[[i]]$abbrevs.y,"group"] <- "foreign"
}
medsLum <- do.call(rbind, medsLum)

ggplot(medsLum) +
  geom_freqpoly(aes(x = median, colour = group), stat = "density") +
  facet_wrap(~qtile,nrow = 4)


# local v nonlocal by sub



