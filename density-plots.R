source("loading-cleaning.r")

df$group <- "both"
df$group[df$substrate == "a"] <- "a"

png(file='plott.png', width=1400, height=1800, units = "px", res = 150)
par(mfrow=c(12,5), mar=c(.3,.3,.3,.3))
for(i in 1:length(unique(df$pop_spp))){
  fltB <- df[df$pop_spp == unique(df$pop_spp)[i] & df$group == "both",]
  fltA <- df[df$pop_spp == unique(df$pop_spp)[i] & df$group == "a",]
  plot(1,xlim = c(0,0.75), ylim = c(0,18),
       xlab = "lumMean", ylab = "density",
       xaxt = 'n', yaxt = 'n')
  title(unique(df$abbrevs)[i], adj = 0.9, line = -1, cex.main = 0.8)
  lines(density(fltB$lumMean), col = "black", lty = 2, lwd = 2)
  lines(density(fltA$lumMean), col = "black")
  lines(density(fltB[fltB$substrate =="c",]$lumMean), col = alpha(rgb(1,0,0), 0.8))
  lines(density(fltB[fltB$substrate =="b",]$lumMean), col = alpha(rgb(0,0,1), 0.8))
}
dev.off() 

# cute function
plotfn= function(u) {
  flt = filter(df, pop_spp == u)
  ggplot(flt, aes(x = lumMean, colour = group)) +
    geom_density() + theme_bw() 
}
pp <- lapply(unique(df$pop_spp),plotfn)