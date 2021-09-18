source("loading-cleaning.R")
# Are ROIs scaled correctly? YES. Although when scale bar mistake, weird results
# Compare in imageJ
grep("aggen", unique(df$pop_spp), value = T)
aggenys <- df[df$pop_spp == "aggeneys_olivacea_CDE",]
b14 <- aggenys[aggenys$mspec == "CDE_b14",]

# size distributions of rocks across pops
b <- df[df$substrate == "b",]

source("loading-cleaning.r")
library("plotrix")

df$group <- "both"
df$group[df$substrate == "a"] <- "a"

# order by species
sppOrder <- df[, c("abbrevs", "pop_spp")]
sppOrder <- separate(sppOrder, "abbrevs", into = c("pop", "spp", "code"), sep = "_")
sppOrder <- sppOrder[order(sppOrder$spp, decreasing = F),]
unique(df$pop_spp)

df$spp_pop <- df$pop_spp
df <- separate(df, spp_pop, into = c("pop", "spp", "code"), sep = "_")
df <- df[order(df$spp, decreasing = F),]
df <- unite(df, col = "spp_pop", c("spp", "pop", "code"), sep = "_")

df$logArea <- log(df$area)
png(file='output//areaa-distributions.png', width=1400, height=1800, units = "px", res = 120)
par(mfrow=c(12,5), mar=c(.3,.3,.3,.3))
for(i in 1:length(unique(df$pop_spp))){
  fltB <- df[df$pop_spp == unique(df$pop_spp)[i] & df$group == "bg",]
  fltA <- df[df$pop_spp == unique(df$pop_spp)[i] & df$group == "lithops",]
  plot(1,xlim = c(0,10), ylim = c(0,2),
       xlab = "lumMean", ylab = "density",
       xaxt = 'n', yaxt = 'n')
  title(unique(df$abbrevs)[i], adj = 0.9, line = -1, cex.main = 0.8)
  lines(density(fltA$logArea), col = "black")
  lines(density(fltB[fltB$substrate =="c",]$logArea), col = alpha(rgb(1,0,0), 0.8))
  lines(density(fltB[fltB$substrate =="b",]$logArea), col = alpha(rgb(0,0,1), 0.8))
}
dev.off() 

# is there a size brightness association? YES
# insert lithops mean size
# insert abbrevs
# get sample sizes - overall and around lithops size

lsdf <- split(df, df$pop_spp)
lsdf <- lapply(lsdf, function(x) {
  glm(x$lumMean~x$logArea)
})
lsdf <- lapply(lsdf, function(x) {
  coef(x)
})
lsdf <- lapply(lsdf, function(x) {
  x[2]
})
x <- names(lsdf)
y <- unlist(lsdf, use.names = F) 
names(y) <- x
y <- sort(y)

df$pop_spp <- as.factor(df$pop_spp)
df$pop_spp <- factor(df$pop_spp, levels = names(y))
unique(df$pop_spp)

png(file='output//size-bright-corrr.png', width=1400, height=1800, units = "px", res = 150)
par(mfrow=c(12,5), mar=c(.3,.3,.3,.3))
for(i in 1:length(levels(df$pop_spp))){
  fltB <- df[df$pop_spp == levels(df$pop_spp)[i] & df$group == "both",]
  fltA <- df[df$pop_spp == levels(df$pop_spp)[i] & df$group == "a",]
  plot(fltB[fltB$substrate =="b",]$logArea, fltB[fltB$substrate =="b",]$lumMean)
  fit <- glm(fltB$lumMean~fltB$logArea)
  co <- coef(fit)
  abline(fit, col="blue", lwd=2)
  title(unique(fltB$abbrevs), adj = 0.9, line = -1, cex.main = 0.8)
  
  
}
dev.off() 

png(file='output//size-bright-cor.png', width=1400, height=1800, units = "px", res = 150)
par(mfrow=c(12,5), mar=c(1,1,1,1))
for(i in 1:length(levels(df$pop_spp))){
  fltB <- df[df$pop_spp == levels(df$pop_spp)[i] & df$group == "both",]
  fltA <- df[df$pop_spp == levels(df$pop_spp)[i] & df$group == "a",]
  plot(fltB[fltB$substrate =="b",]$area, fltB[fltB$substrate =="b",]$lumMean)
  points(mean(fltA$area,), mean(fltA$lumMean), pch = 21, colour = "red", bg = "red")
  fit <- glm(fltB$lumMean~fltB$area)
  co <- coef(fit)
  abline(fit, col="blue", lwd=2)
  title(unique(fltB$abbrevs), adj = 0.9, line = -1, cex.main = 0.8)
  
  
}
dev.off() 

bavi <- df[df$pop_spp == grep("Baviaans", unique(df$pop_spp), value =  T),]
mean(bavi[bavi$substrate == "a", "area"])
