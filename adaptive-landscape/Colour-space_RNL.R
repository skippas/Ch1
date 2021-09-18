source("loading-cleaning.R")
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)

df <- as.data.table(df)
wf <- c(0.5656854 ,0.1414214 ,0.1000000) # these are your human e values for each cone

# x coord (R:G Opp mech)
df[,lwLog := log(lwMean)]
df[,mwLog := log(mwMean)]
df[,swLog := log(swMean)]

opp1NoiseTerm <- sqrt(1/(wf[2]^2 + wf[3]^2))
df[,xCoord := (lwLog - mwLog)*opp1NoiseTerm]

# y coord 
wfC <- combn(wf,2)
opp2NoiseTerm <- sqrt((wf[2]^2 + wf[3]^2) / 
                        (prod(wfC[,1])^2 + prod(wfC[,2])^2 + prod(wfC[,3])^2) )
df[,yCoord := (df[,lwLog*(wf[2]^2 / (wf[2]^2 + wf[3]^2))] + 
                        df[,mwLog*(wf[3]^2 / (wf[2]^2 + wf[3]^2))])]
df[,yCoord := (swLog*opp2NoiseTerm) - (yCoord*opp2NoiseTerm)]
#df[,yCoord := yCoord*opp2NoiseTerm]


# remove colour outliers
plot(df$xCoord, df$yCoord, abline(h=c(0.5,-3.5), v = 4))
# conditions
outliers <- df[df$yCoord < -3.5 | df$yCoord > 0.5 | df$xCoord > 4]
df <- df[df$yCoord > -3.5 & df$yCoord < 0.5 & df$xCoord < 4]
plot(df$xCoord, df$yCoord, abline(h=c(0.5,-3.5), v = 4))

# plots in RNL space
# order by species
df$spp_pop <- df$pop_spp
df <- separate(df, spp_pop, into = c("pop", "spp", "code"), sep = "_")
df <- df[order(df$spp, decreasing = F),]
df <- unite(df, col = "spp_pop", c("spp", "pop", "code"), sep = "_")

# ggplots
df <- df[df$substrate == "a" | df$substrate == "b" | df$substrate == "c",]
plot_lst <- vector("list", length = length(unique(df$pop_spp)))

# Set a colour palette
library(RColorBrewer)
myColors <- brewer.pal(3,"Set1")
names(myColors) <- c("c", "b", "a")
colScale <- scale_colour_manual(name = "substrate",values = myColors)
colScaleFill <- scale_fill_manual(name = "substrate",values = myColors)

for(i in 1:length(unique(df$pop_spp))){ # claus wilke solution
   
    xx <- df[df$pop_spp == unique(df$pop_spp)[i],]
  
  pmain <- ggplot(xx, aes(x = xCoord, y = yCoord, color = substrate)) +
    geom_point(size = 0.5, stroke = 0, shape = 16) +  xlim(-1,4) + ylim(-3.5,0.5) +
    geom_point(x=0, y=0, shape = 3, color = "black") +
    theme(axis.title=element_blank(),
          #axis.text=element_blank(),
          #axis.ticks=element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = "none") +
    annotate("text", x = 2, y = 0.5,size = 2, label = unique(xx$abbrevs)) +
    annotate("text", x = 1, y = -3.5,size = 2, label = unique(xx$pop_spp)) +
    colScale
  # Marginal densities along x axis
  xdens <- axis_canvas(pmain, axis = "x")+
    geom_density(data = xx, aes(x = xCoord, fill = substrate),
                 alpha = 0.7, size = 0.2) + 
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    colScaleFill
  # Marginal densities along y axis
  # Need to set coord_flip = TRUE, if you plan to use coord_flip()
  ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
    geom_density(data = xx, aes(x = yCoord, fill = substrate),
                 alpha = 0.7, size = 0.2)+
    coord_flip()+
    colScaleFill
  p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  ggdraw(p2)
  
  plot_lst[[i]] <- p2
}

plotsPerPng <- 1:length(plot_lst)
plotsPerPng <- split(plotsPerPng, ceiling(seq_along(plotsPerPng)/9))
for(i in 1:length(plotsPerPng)){
  png(file= paste0('output//adaptive-landscape//human-RNL-space', i, ".png"), width=1400, height=1800, units = "px", res = 300)
  p <- cowplot::plot_grid(plotlist = plot_lst[plotsPerPng[[i]]], nrow = 3, ncol = 3)
  print(p)
  dev.off() 
}

# Geometric mean positions
geomMeans <- df %>% group_by(abbrevs, substrate) %>%
  summarise_at(vars(swLog, mwLog, lwLog), mean)

geomMeans <- as.data.table(geomMeans)
# x coord (R:G Opp mech)
opp1NoiseTerm <- sqrt(1/(wf[2]^2 + wf[3]^2))
geomMeans[,xCoord := (lwLog - mwLog)*opp1NoiseTerm]

# y coord 
wfC <- combn(wf,2)
opp2NoiseTerm <- sqrt((wf[2]^2 + wf[3]^2) / 
                        (prod(wfC[,1])^2 + prod(wfC[,2])^2 + prod(wfC[,3])^2) )
geomMeans[,yCoord := (geomMeans[,lwLog*(wf[2]^2 / (wf[2]^2 + wf[3]^2))] + 
                        geomMeans[,mwLog*(wf[3]^2 / (wf[2]^2 + wf[3]^2))])]
geomMeans[,yCoord := (swLog*opp2NoiseTerm) - (yCoord*opp2NoiseTerm)]

png(file= paste0('output//adaptive-landscape//geom-means_human-RNL-space.png'), 
    width=1400, height=1800, units = "px", res = 400)
ggplot(geomMeans, aes(x = xCoord, y = yCoord, color = substrate)) +
  geom_point(size = 1, stroke = 0.5, shape = 1) +  xlim(-0.2,3) + ylim(-2.5,0.2) +
  geom_point(x=0, y=0, shape = 3, color = "black", stroke = 1) +
  theme(axis.title=element_blank(),
        #axis.text=element_blank(),
        #axis.ticks=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.position = "none") + theme_bw()+
  colScale
dev.off()

# geologies
geol <- df[df$substrate == "b",]
geol <- geol[geol$geology != "graniteOrComplex"]

p <- ggplot(geol, aes(x = xCoord, y = yCoord, color = geology)) +
  geom_point(size = 1, stroke = 0.5, shape = 1) +  
  geom_point(x=0, y=0, shape = 3, color = "black", stroke = 1) +
  stat_ellipse(color = "black") +
  theme(axis.title=element_blank(),
        #axis.text=element_blank(),
        #axis.ticks=element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        legend.position = "none") + theme_bw()

png(file= paste0('output//adaptive-landscape//geol-chrom.png'), 
    width=1400, height=1800, units = "px", res = 400)
p + facet_grid(rows = vars(geology))
dev.off()

# BRIGHTNESS




xx <- df[df$pop_spp == unique(df$pop_spp)[3],]
plot(xx$xCoord, xx$yCoord)

plot(1, xlim = c(0,15), ylim = c(-20,0))
points(xx[xx$substrate == "a"]$xCoord,xx[xx$substrate == "a"]$yCoord, pch = 20, col = alpha(rgb(.3,.3,.3),0.8))
points(xx[xx$substrate == "b"]$xCoord,xx[xx$substrate == "b"]$yCoord, pch = 16, col = alpha(rgb(0,0,1), 0.2))
points(xx[xx$substrate == "c"]$xCoord,xx[xx$substrate == "c"]$yCoord, pch = 16, col = alpha(rgb(1,0,0), 0.2))


# cross reffing coords with pavo colour dists to check they make sense

x <- df[df$pop_spp == grep("kenh", unique(df$pop_spp), value = T),]
x <- x[x$mspec == grep("1668", unique(x$mspec), value = T),]
x <- x[x$roi == "b6",]
naam <- paste0(x$abbrevs,x$mspec,x$roi)
x <- x[, c("swMean", "mwMean", "lwMean","lumMean")]
row.names(x) <- naam

y <- df[df$pop_spp == grep("Haasriver2", unique(df$pop_spp), value = T),]
y <- y[y$mspec == grep("3788", unique(y$mspec), value = T),]
y <- y[y$roi == "b6",]
naam <- paste0(y$abbrevs,y$mspec,y$roi)
y <- y[, c("swMean", "mwMean", "lwMean","lumMean")]
row.names(y) <- naam

z <- df[df$pop_spp == grep("boeg", unique(df$pop_spp), value = T),]
z <- z[z$mspec == grep("1952", unique(z$mspec), value = T),]
z <- z[z$roi == "a2",]
naam <- paste0(z$abbrevs,z$mspec,z$roi)
z <- z[, c("swMean", "mwMean", "lwMean","lumMean")]
row.names(y) <- naam

yz <- rbind(y,z)
colnames(yz) <- c("s", "m", "l", "lum") ; x
coldist(yz, n = c(1,16,32), weber = .1, weber.achro = .1, 
         qcatch = "Qi", achromatic = TRUE)

data("sicalis")

vismod1 <- vismodel(sicalis,
                    visual = "canis", scale = 10000,
                    illum = "D65", relative = FALSE
)




