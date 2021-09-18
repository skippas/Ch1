source("loading-cleaning.R")
library(data.table)
library(patchwork)
library(ggpmisc)

# sample single Lithops 

# Theme set
theme_set(theme_bw())

library(RColorBrewer)
#myColors <- brewer.pal(3,"Set1")
myColors <- c('#696969', "#00BFC4", "#F8766D")
names(myColors) <- c("a", "b", "c")
colScale <- scale_colour_manual(name = "substrate",
                                values = myColors,
                                labels = c('Lithops','Rock', 'Soil'),
                                drop = T)
colScaleFill <- scale_fill_manual(name = "substrate",
                                  values = myColors,
                                  labels = c('Lithops','Rock', 'Soil'),
                                  drop = T)

# substrates means in visual space across populations ####

source("functions//gmean.R")
# Luminance geometric means
subGmeans<- df %>% group_by(abbrevs, substrate) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) 

meanLumP <- ggplot(subGmeans, aes(lumMean)) +
  geom_histogram(aes(fill = substrate), color = "black") + 
  theme_minimal()+ facet_wrap(vars(substrate), dir = "v") +
  theme(strip.text = element_blank(), legend.position = 'none') +
  xlab(label = "Geometric mean luminance")+
  scale_y_continuous(labels = scales::label_number(accuracy = 1))+
  colScaleFill

# Colour
source("functions//RNL_colspace_tri.R")
subGmeans <- RNL_colspace_tri(df = subGmeans)

meanChromaP <-ggplot(subGmeans, aes(x = xCoord, y = yCoord,
                      colour = substrate, shape = substrate))+
  geom_point()+
  xlab("G:R (JNDs)")+
  ylab("Y:B (JNDs)")+
  colScale+
  scale_shape_manual(values = c(16,15,17),
                     labels = c('Lithops', 'Rock', 'Soil'))+
  theme(legend.position = 'top')

png("output//adaptive-landscape//RvS//Gmean-col-scatter-bySub.png", 
      type = 'cairo',
      width=200, height=100, units = "mm", res = 300)
meanLumP + meanChromaP + plot_layout(ncol = 2)
dev.off()


rm(subGmeans, meanLumP, meanChromaP)

# local substrate mean distances ####

###* distribution of local distances betw R and S ####
# doesn't tell the whole story because of spread
# TO DO: local dL distances when nearest subsets used
source("functions//gMeanDists.R")
gMeanDistDF <- gMeanDists(df = df)
locDist <- gMeanDistDF[gMeanDistDF$comparison == "local",]
locDistbc <- locDist[locDist$patch1 %in% c("b", "c") & patch2 %in% c("b", "c"),]
locDistbc <- pivot_longer(locDistbc, cols = c("dL", "dS"),
                          names_to = "visfeature", values_to =  "gMeanDist")

png("output//adaptive-landscape//gMeanDists-betw-locSubs.png", 
    width=800, height=600, units = "px", res = 100)
ggplot(locDistbc, aes(gMeanDist)) +
  geom_histogram(aes(fill = visfeature)) + 
  facet_wrap(vars(visfeature), dir = "v", scales = "free")+ 
  theme(strip.text = element_blank())+
  theme_minimal()+
  xlab("Local rock-soil geometric mean contrast (JNDs)")+ 
  scale_fill_manual(name = "Visual feature",values = c("dL" = "black", "dS" = "coral2"),
                     labels = c("Luminance", "Chroma"))
dev.off()

###* distances between local substrates per pop ####
# library(urbnthemes)
# set_urbn_defaults(style = 'print')

visInfoRename <- c(dL = 'Lum dist', dS = 'Chroma dist')
png("output\\adaptive-landscape\\RvS\\RS-distances-perPop-LumCol.png",
    type = 'cairo', units = 'mm',
    height = 300, width = 200, res = 300)
locDistbc %>% ungroup() %>%
  mutate(visfeature = as.factor(visfeature),
         bgpop = reorder_within(bgpop, gMeanDist, visfeature)) %>%
  ggplot(aes(bgpop, gMeanDist))+
  geom_col(show.legend = FALSE)+
  facet_wrap(~visfeature, scales = "free_y", ncol = 1,
             labeller = labeller(visfeature = visInfoRename))+
  coord_flip()+
  scale_x_reordered()+
  scale_y_continuous(expand = c(0,0))
dev.off()

# Lollipop charts a mission to facet, but possible. see: https://stackoverflow.com/questions/52214071/how-to-order-data-by-value-within-ggplot-facets
png("output\\adaptive-landscape\\RvS\\RS-spread-perPop-colour.png",
    type = 'cairo', units = 'mm',
    height = 300, width = 200, res = 300)
locDistbc %>% filter(visfeature == "dL") %>%
  mutate(bgpop = fct_reorder(bgpop, gMeanDist)) %>%
  ggplot(aes(gMeanDist, bgpop))+
  geom_segment(aes(x = 0, xend = gMeanDist, y = bgpop, yend = bgpop))+ 
  geom_point()+
  scale_x_continuous(expand = expand_scale(mult = c(0, 0)), limits = c(0, 15))+
  labs(x = NULL, y = "Miles Per Gallon")+ theme_urbn_print()

png("output\\adaptive-landscape\\RvS\\RS-spread-perPop-colour.png",
    type = 'cairo', units = 'mm',
    height = 300, width = 200, res = 300)
locDistbc %>% filter(visfeature == "dS") %>%
  mutate(bgpop = fct_reorder(bgpop, gMeanDist)) %>%
  ggplot(aes(gMeanDist, bgpop)) +
  geom_segment(aes(x = 0, xend = gMeanDist, y = bgpop, yend = bgpop)) + 
  geom_point() +
  scale_x_continuous(expand = expand_scale(mult = c(0, 0)), limits = c(0, 2)) +
  labs(x = NULL, y = "Miles Per Gallon")+ theme_urbn_print()

rm(locDist)

# local substrate correlation ####
# urbnthemes::undo_urbn_defaults()
statPoly <-  stat_poly_eq(formula = y~x,
                          aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                            sep = "~~~")),
                          parse = TRUE, vstep = 3, size =2)

source("rock-v-soil//rock-estimate.R")
props <- df %>% count(abbrevs, substrate)
props <- props %>% pivot_wider(names_from = substrate, values_from = n)
props$rock_cover <- (props$b / (props$b + props$c)*100)
props$soil_cover <- (props$c / (props$b + props$c)*100)
props <- arrange(props, as.factor(abbrevs))
source("loading-cleaning.R")

meanCoords <- df %>%
  group_by(abbrevs, substrate) %>%
  summarise_at(vars(xCoord, yCoord, lumMean), mean) %>%
#  pivot_longer(cols = c(xCoord, yCoord, lumMean), # To restore the facet way
#               names_to = "visInfo") %>% 
  pivot_wider(names_from = substrate,
              values_from = c(xCoord, yCoord,lumMean)) %>%
#  select(-c(a)) %>%
  cbind(., props[match(.$abbrevs, props$abbrevs),
                 c("rock_cover", "soil_cover")]) 

p1<- ggplot(meanCoords,
            aes(lumMean_b, lumMean_c, colour = rock_cover))+
  geom_point()+
  geom_abline(lty=2)+
#  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .5)+
  coord_fixed(ratio = 1)+
  scale_x_continuous(expand = c(0,0), limits = c(0,.5))+
  scale_y_continuous(expand = c(0,0), limits = c(0,.5))+
  xlab("R mean lum")+
  ylab("S mean lum")+
  statPoly+
  theme(legend.position = 'none')

p2<- ggplot(meanCoords,
            aes(xCoord_b, xCoord_c, colour = rock_cover))+
  geom_point()+
  geom_abline(lty=2)+
#  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .5)+
  coord_fixed(ratio = 1)+
  scale_x_continuous(expand = c(0,0), limits = c(0,3))+
  scale_y_continuous(expand = c(0,0), limits = c(0,3))+
  xlab("R mean G:R")+
  ylab("S mean G:R")+
  statPoly+
  theme(legend.position = 'none')

p3<- ggplot(meanCoords,
            aes(yCoord_b, yCoord_c, colour = rock_cover))+
  geom_point()+
  geom_abline(lty=2)+
#  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .5)+
  coord_fixed(ratio = 1)+
  scale_x_continuous(expand = c(0,0), limits = c(-2.2,-.3))+
  scale_y_continuous(expand = c(0,0), limits = c(-2.2,-.3))+
  xlab("R mean Y:B")+
  ylab("S mean Y:B")+
  labs(colour = 'Rock percent')+
  statPoly
  #theme(legend.position = 'none')

lsLumChromRSMeanRelP <- list(p1,p2,p3)


# dummy <- df %>%
#   group_by(abbrevs, substrate) %>%
#   summarise_at(vars(xCoord, yCoord, lumMean), mean) 
# 
# dummy <- data.frame(b = c(range(dummy$lumMean),
#                           range(dummy$xCoord),
#                           range(dummy$yCoord)),
#                     c = c(range(dummy$lumMean),
#                           range(dummy$xCoord),
#                           range(dummy$yCoord)),
#                     visInfo = c(rep("lumMean",2),
#                                 rep("xCoord",2),
#                                 rep("yCoord",2)),
#                     rock_cover = 0)
# 
# facetLab = c(lumMean = 'Lum', xCoord = 'RG', yCoord = 'BY')
# library(ggpmisc)
# 
# lumChromRSMeanRelP <- 
#   ggplot(meanCoords, aes(b,c, colour = rock_cover))+
#  geom_point()+
#   geom_blank(data = dummy)+
#   geom_abline(lty = 2)+
#   coord_cartesian()+
#   theme_minimal()+
#   theme(aspect.ratio = 1,
#         legend.position = 'top')+
#   xlab("Rock mean")+
#   ylab("Soil mean")+
#   stat_poly_eq(formula = y~x,
#                aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
#                                  sep = "~~~")),
#                parse = TRUE, vstep = 3, size =3)+
#   facet_wrap(vars(visInfo),scales = 'free', nrow = 1,
#              labeller = labeller(visInfo = facetLab))
# png("output\\adaptive-landscape\\RvS\\RS-correlation.png",
#     type = 'cairo', units = 'mm',
#     height = 150, width = 215, res = 300)
# lumChromRSMeanRelP
# dev.off()

# Substrate heterogeneity ####
###* MAD ####
###** Plot dataframes ####
source("functions//MAD.R")
source('functions//RNL_colspace_tri.R')
  
  
colSpread <- RNL_colspace_tri(df = df)
colSpread <- colSpread %>% group_by(abbrevs, substrate) %>% 
  summarize(MADXcoord = MAD(xCoord), MADYcoord = MAD(yCoord))
colSpread <- pivot_wider(colSpread, names_from = "substrate",
                         values_from = c("MADXcoord", "MADYcoord")) 
colSpreadDist <- colSpread
colSpreadDist$gMeanDist <- locDistbc[match(colSpread$abbrevs,
                                           locDistbc$bgpop),]$gMeanDist

lumSpread <- df %>% group_by(abbrevs, substrate) %>% 
  summarize(MADlum = MAD(lumMean))
lumSpread <- pivot_wider(lumSpread, names_from = "substrate",
                         values_from = c("MADlum")) 
lumSpreadDist <- merge(lumSpread, locDistbc[locDistbc$visfeature == "dL",c("gMeanDist", "lithpop")], 
                       by.x = "abbrevs", by.y = "lithpop")
textLabs <- lumSpreadDist[lumSpreadDist$c / lumSpreadDist$b > 2 |
                            lumSpreadDist$c / lumSpreadDist$b < 0.5,]

###** Plots ####
###*** Scatterplots relship RS spread ####
library(ggrepel)

p1<-ggplot(lumSpreadDist, aes(b, c, colour = gMeanDist))+
  geom_point()+
  coord_fixed(ratio = 1)+
  geom_abline(lty=2)+
#  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .5)+
  scale_x_continuous(expand = c(0,0), limits = c(0,0.16))+
  scale_y_continuous(expand = c(0,0), limits = c(0,.16))+
  scale_color_viridis_c(option = 'magma')+
  xlab("R MAD lum")+
  ylab("S MAD lum")+
  statPoly+
  theme(legend.position = 'none')

p2<- ggplot(colSpreadDist, aes(MADXcoord_b, MADXcoord_c, 
                          colour = gMeanDist, label = abbrevs))+
  geom_point()+
  geom_abline(lty=2)+
  coord_fixed(ratio = 1, clip = 'off')+
#  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .5)+
  scale_x_continuous(expand = c(0,0), limits = c(0,.75))+
  scale_y_continuous(expand = c(0,0), limits = c(0,.75))+
  scale_color_viridis_c(option = 'magma')+
  xlab("R MAD G:R")+ 
  ylab("S MAD G:R")+
  statPoly+
  theme(legend.position = 'none')

p3<-ggplot(colSpreadDist, aes(MADYcoord_b, MADYcoord_c, 
                          colour = gMeanDist, label = abbrevs))+
  geom_point()+
  geom_abline(lty=2)+
#  stat_summary(fun.data= mean_cl_normal)+ 
  geom_smooth(method='lm', colour = 'grey70', se = F, lwd = .5)+
  coord_fixed(ratio = 1)+
  scale_x_continuous(expand = c(0,0), limits = c(0,0.6))+
  scale_y_continuous(expand = c(0,0), limits = c(0,.6))+
  scale_color_viridis_c(option = 'magma')+
  xlab("R MAD Y:B ")+
  ylab("S MAD Y:B")+
  statPoly+
  labs(colour = 'R-S mean dist')

png("output\\adaptive-landscape\\RvS\\RS-spread-relship-scatterp.png",
    type = 'cairo', units = 'mm',
    height = 290, width = 215, res = 300)
p1 + p2 + p3 + plot_layout(ncol = 3)
dev.off()

png("output\\adaptive-landscape\\RvS\\RS-spread&mean-relship-scatterp.png",
    type = 'cairo', units = 'mm',
    height = 140, width = 215, res = 300)
lsLumChromRSMeanRelP[[1]] + lsLumChromRSMeanRelP[[2]] +
  lsLumChromRSMeanRelP[[3]] + p1 + p2 + p3 +plot_layout(ncol = 3)
dev.off()

rm(lsLumChromRSMeanRelP, statPoly)

###*** Bar charts RS spread per population ####
library(tidytext) # https://juliasilge.com/blog/reorder-within/

colSpread <- colSpread[,c("MADXcoord_b","MADXcoord_c",
                          "MADYcoord_b","MADYcoord_c" , "abbrevs")] %>%
  pivot_longer(cols = c(contains("_b"), contains("_c")),
               names_to = 'visInfo', values_to = 'MAD') 
  #  separate(visInfo, into = c("visInfo", "substrate"), sep="_(?=[^_]+$)")  %>%
  #colSpread$visInfo <- gsub("MAD", "",colSpread$visInfo) 

png("output\\adaptive-landscape\\RvS\\RS-spread-perPop-colour.png",
    type = 'cairo', units = 'mm',
    height = 300, width = 200, res = 300)
colSpread %>% ungroup() %>%
  mutate(visInfo = as.factor(visInfo),
         abbrevs = reorder_within(abbrevs, MAD, visInfo)) %>%
  ggplot(aes(abbrevs, MAD))+
  geom_col(show.legend = FALSE)+
  facet_wrap(~visInfo, scales = "free_y")+
  coord_flip()+
  scale_x_reordered()+
  scale_y_continuous(expand = c(0,0))
dev.off()

lumSpread <- lumSpread[,c("b","c", "abbrevs")] %>%
  pivot_longer(cols = c("c", "b"), names_to = 'visInfo', values_to = 'MAD') 

png("output\\adaptive-landscape\\RvS\\RS-spread-perPop-Lum.png",
    type = 'cairo', units = 'mm',
    height = 300, width = 200, res = 300)
ilumSpread %>% ungroup() %>%
  mutate(visInfo = as.factor(visInfo),
         abbrevs = reorder_within(abbrevs, MAD, visInfo)) %>%
  ggplot(aes(abbrevs, MAD)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~visInfo, scales = "free_y", ncol = 1) +
  coord_flip() +
  scale_x_reordered() +
  scale_y_continuous(expand = c(0,0))
dev.off()



###* IQR of luminance and colour for populations ####
###* difference in heterogeneity between local substrates by population ####



# Distributions of all raw data in col / luminance space ####

###* facet scatterplots of L, R and L, S in colour space ####

theme_set(theme_minimal())
# Reorder by collapseSpp
df <- df %>% 
  ungroup() %>%
  arrange(collapseSpp, abbrevs) %>%
  mutate(abbrevs = factor(abbrevs, levels = unique(abbrevs)))

# LS scatter
png("output\\adaptive-landscape\\RvS\\LS-chroma-scatter-facet-1.png",
    type = 'cairo', units = 'mm', width = 215, height = 300, res = 300)
df[df$abbrevs %in% c(levels(df$abbrevs)[1:28]),] %>% 
  filter(substrate != "b") %>% 
  ggplot(aes(xCoord,yCoord,fill = substrate, colour = substrate))+
  geom_point(pch = 21, size = 1, alpha =.3)+
  facet_wrap(vars(abbrevs), scales = 'free', ncol = 4)+
  theme(legend.position = 'top')
dev.off()
png("output\\adaptive-landscape\\RvS\\LS-chroma-scatter-facet-2.png",
    type = 'cairo', units = 'mm', width = 215, height = 280, res = 300)
df[df$abbrevs %in% c(levels(df$abbrevs)[29:56]),] %>% 
  filter(substrate != "b") %>% 
  ggplot(aes(xCoord,yCoord,fill = substrate, colour = substrate))+
  colScale+ colScaleFill+
  geom_point(pch = 21, size = 1, alpha =.3)+
  facet_wrap(vars(abbrevs), scales = 'free', ncol = 4)+
  theme(legend.position = 'none')
dev.off()

# LR Scatter
png("output\\adaptive-landscape\\RvS\\LR-chroma-scatter-facet-1.png",
    type = 'cairo', units = 'mm', width = 200, height = 300, res = 300)
df[df$abbrevs %in% c(levels(df$abbrevs)[1:28]),] %>% 
  filter(substrate != "c") %>% 
  ggplot(aes(xCoord,yCoord,fill = substrate, colour = substrate))+
  geom_point(pch = 21, size = 1, alpha =.3)+
  facet_wrap(vars(abbrevs), scales = 'free', ncol = 4)
dev.off()
png("output\\adaptive-landscape\\RvS\\LR-chroma-scatter-facet-2.png",
    type = 'cairo', units = 'mm', width = 200, height = 300, res = 300)
df[df$abbrevs %in% c(levels(df$abbrevs)[29:56]),] %>% 
  filter(substrate != "c") %>% 
  ggplot(aes(xCoord,yCoord,fill = substrate, colour = substrate))+
  geom_point(pch = 21, size = 1, alpha =.3)+
  facet_wrap(vars(abbrevs), scales = 'free', ncol = 4)
dev.off()

# RS scatter
png("output\\adaptive-landscape\\RvS\\RS-chroma-scatter-facet-1.png",
    type = 'cairo', units = 'mm', width = 215, height = 300, res = 300)
df[df$abbrevs %in% c(levels(df$abbrevs)[1:28]),] %>% 
  filter(substrate != "a") %>% 
  ggplot(aes(xCoord,yCoord,fill = substrate, colour = substrate))+
  geom_point(pch = 21, size = 1, alpha =.3)+
  facet_wrap(vars(abbrevs), scales = 'free', ncol = 4)+
  scale_fill_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                    labels = c("Rock", "Soil"))+
  scale_colour_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                      labels = c("Rock", "Soil"))+
  ylab('Yellow:Blue opponent mechanism')+
  xlab('Green:Red opponent mechanism')+
  theme(legend.position = 'top')
dev.off()
png("output\\adaptive-landscape\\RvS\\RS-chroma-scatter-facet-2.png",
    type = 'cairo', units = 'mm', width = 215, height = 280, res = 300)
df[df$abbrevs %in% c(levels(df$abbrevs)[29:56]),] %>% 
  filter(substrate != "a") %>% 
  ggplot(aes(xCoord,yCoord,fill = substrate, colour = substrate))+
  scale_fill_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                    labels = c("Rock", "Soil"))+
  scale_colour_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                    labels = c("Rock", "Soil"))+
  ylab('Yellow:Blue opponent mechanism')+
  xlab('Green:Red opponent mechanism')+
  geom_point(pch = 21, size = 1, alpha =.3)+
  facet_wrap(vars(abbrevs), scales = 'free', ncol = 4)+
  theme(legend.position = 'none')
dev.off()

###* facet plots in colour space WITH marginal densities ####
source("functions//RNL_colspace_tri.R")
df <- RNL_colspace_tri(df = df)

# remove colour outliers
plot(df$xCoord, df$yCoord, abline(h=c(0.5,-3.5), v = 4))
# conditions
outliers <- df[df$yCoord < -3.5 | df$yCoord > 0.5 | df$xCoord > 4]
df <- df[df$yCoord > -3.5 & df$yCoord < 0.5 & df$xCoord < 4]
plot(df$xCoord, df$yCoord, abline(h=c(0.5,-3.5), v = 4))

# plots in RNL space
# order by species
df$species <- fct_relevel(df$species, sort)

# ggplots
plot_lst <- vector("list", length = length(unique(df$pop_spp)))
# Set a colour palette
library(RColorBrewer)
myColors <- brewer.pal(3,"Set1")
names(myColors) <- c("c", "b", "a")
colScale <- scale_colour_manual(name = "substrate",values = myColors)
colScaleFill <- scale_fill_manual(name = "substrate",values = myColors)

library(cowplot)
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
  xdens <- axis_canvas(pmain, axis = "x") +
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
  png(file= paste0('output//adaptive-landscape//human-RNL-space', i, ".png"),
      width=1400, height=1800, units = "px", res = 300)
  p <- cowplot::plot_grid(plotlist = plot_lst[plotsPerPng[[i]]], nrow = 3, ncol = 3)
  print(p)
  dev.off() 
}
pdf("output//adaptive-landscape//human-RNL-colspace.pdf", onefile = T)
for(i in 1:length(plotsPerPng)){
  p <- cowplot::plot_grid(plotlist = plot_lst[plotsPerPng[[i]]], nrow = 3, ncol = 3)
  print(p)
}
dev.off() 


###* facet plots in luminance space ####

library("plotrix")

df$logLum <- exp(log(df$lumMean)) # whats the math effect of doing this!!!
sumStats <- df %>% group_by(pop_spp) %>% filter(substrate == "b") %>%
  summarise(rock_sd = sd(logLum), rock_n = n())

df$abbrevs <- fct_relevel(df$abbrevs, levels = unique(df[order(df$species), "abbrevs"]))

par(mfrow=c(12,5), mar=c(.3,.3,.3,.3))
for(i in 1:length(levels(df$abbrevs))){
  plot(-10,xlim = c(0,1), ylim = c(0,18),
       xlab = "logLum", ylab = "density",
       xaxt = 'n', yaxt = 'n')
  title(levels(df$abbrevs)[i], adj = 0.9, line = -1, cex.main = 0.8)
  title(unique(df[df$abbrevs == levels(df$abbrevs)[i],"geology"]), 
        adj = 0.9, line = -2, cex.main = 0.8)
  lines(density(df[df$abbrevs == levels(df$abbrevs)[i] & df$substrate == "b",]$logLum), 
        col = "blue", lty = 1, lwd = 1)
  lines(density(df[df$abbrevs == levels(df$abbrevs)[i] & df$substrate == "c",]$logLum), 
        col = "red", lty = 1, lwd = 1)
}

logLumPanels <- recordPlot()
png(file='output//adaptive-landscape//loglum-dens_by-site.png', width=1400,
    height=1800, units = "px", res = 150)
logLumPanels
dev.off() 
pdf(file='output//adaptive-landscape//loglum-dens_by-site.pdf')
logLumPanels
dev.off() # check out pryr for saving code as object and repeating it

###* facet distributions each pop R S luminance violin plots ####
#function from
#https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
library(gtable)
library(cowplot)

shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}

p <- ggplot(df %>% filter(substrate == c("b", "c")), aes(x= abbrevs, y = lumMean)) + 
  geom_violin(aes(fill = substrate))
png(file= paste0('output//adaptive-landscape//rock-soil-lum-distrib.png'), 
    width=200, height=300, units = "mm", res = 300, type = 'cairo')
p <- p + facet_wrap(facets = vars(abbrevs), scales = "free", ncol = 6) + 
  theme(strip.background = element_blank(), strip.text = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0, "lines")) +
  scale_fill_viridis_d(begin = 0,end = 0.4, direction = -1, labels = c("Rock", "Soil"))
p <- p +
  guides(fill = guide_legend(title.position = "left",
                             label.position = "bottom",
                             title.vjust = 0.7,
                             nrow = 1)) + 
  theme(legend.direction = "horizontal")
grid::grid.draw(shift_legend(p))
dev.off()
# TO DO: ridgeline plots of the diff species, their variation across pops




###* facet luminance histograms R S ####

plotF <- function(start = 1, end = length(unique(df$abbrevs))){
  df %>% 
    filter(as.integer(abbrevs) %in% start:end & substrate != "a") %>%
  ggplot(aes(x=lumMean))+ 
    stat_bin(aes(fill = substrate),
             position="identity",
             alpha=0.5)+
    facet_wrap(~abbrevs, scales = 'free_y',ncol = 4)+
    theme_bw()
}

png("output\\adaptive-landscape\\RvS\\facet-lum-hist-RS_1.png",
    type = 'cairo', units = 'mm', width = 200, height = 300, res = 300)
print(plotF(start = 1, end = 28))
dev.off()
png("output\\adaptive-landscape\\RvS\\facet-lum-hist-RS_2.png",
    type = 'cairo', units = 'mm', width = 200, height = 300, res = 300)
print(plotF(start = 29, end = length(unique(df$abbrevs))))
dev.off()

png("output\\adaptive-landscape\\RvS\\facet-lum-hist-RS_All.png",
    type = 'cairo', units = 'mm', width = 200, height = 280, res = 300)
df %>% 
  filter(as.integer(abbrevs) %in% 1:length(unique(df$abbrevs)) &
           substrate != "a") %>%
  ggplot(aes(x=lumMean))+ 
  stat_bin(aes(fill = substrate),
           position="identity",
           alpha=0.5)+
  scale_fill_manual(values = c("b" = "#00BFC4", "c" = "#F8766D"),
                    labels = c("Rock", "Soil"))+
  xlab('Luminance')+
  facet_wrap(~abbrevs, scales = 'free_y',ncol = 5)+
  theme_bw()+
  theme(legend.position = 'top')
dev.off()
###* Facet density plots ####

library(scales)
show_col(hue_pal()(2))
myColors <- c('grey10', '#F8766D', '#00BFC4')
names(myColors) = c('Lithops','Rock', 'Soil')
colScale <- scale_color_manual(name = 'subRename', values = myColors)

png("output\\adaptive-landscape\\RvS\\facet-lum-density.png",
    type = 'cairo', units = 'mm', width = 215, height = 300, res = 300)
df %>%
  filter(as.integer(abbrevs) %in%
           1:length(unique(df$abbrevs))) %>%
  ggplot(aes(x=lumMean, colour = subRename))+ 
  geom_density()+
  colScale+
  facet_wrap(~abbrevs, scales = 'free_y',ncol = 5)+
  theme_bw()+
  #theme(legend.position = 'none')+
  scale_x_continuous(limits = c(0,0.7))+
  xlab('Luminance')
dev.off()

png("output\\adaptive-landscape\\RvS\\facet-RG-density.png",
    type = 'cairo', units = 'mm', width = 215, height = 300, res = 300)
df %>% 
  filter(as.integer(abbrevs) %in% 1:length(unique(df$abbrevs))) %>%
  ggplot(aes(x=xCoord, colour = subRename))+ 
  geom_density()+
  colScale+
  facet_wrap(~abbrevs, scales = 'free_y',ncol = 5)+
  theme_bw()+
  #theme(legend.position = 'none')+
  #scale_x_continuous(limits = c(0,0.7))+
  xlab('R:G opponent')
dev.off()

png("output\\adaptive-landscape\\RvS\\facet-BY-density.png",
    type = 'cairo', units = 'mm', width = 215, height = 300, res = 300)
df %>% 
  filter(as.integer(abbrevs) %in% 1:length(unique(df$abbrevs))) %>%
  ggplot(aes(x=yCoord, colour = subRename))+ 
  geom_density()+
  colScale+
  facet_wrap(~abbrevs, scales = 'free_y',ncol = 5)+
  theme_bw()+
  #theme(legend.position = 'none')+
  #scale_x_continuous(limits = c(0,0.7))+
  xlab('B:Y opponent')
dev.off()

x<- df %>% 
  filter(subRename != 'Lithops') %>%
  add_count(abbrevs, name = 'n') %>%
  add_count(abbrevs, subRename, name = 'nSubs') %>%
  mutate(subProp = nSubs / n)
df<- merge(df, unique(x[,c('abbrevs','substrate','subProp')]),
           all.x = T) 
df$subProp[is.na(df$subProp)] <- 1

png("output\\adaptive-landscape\\RvS\\facet-lum-density-weighted.png",
    type = 'cairo', units = 'mm', width = 215, height = 300, res = 300)
df %>%
  filter(as.integer(abbrevs) %in%
           1:length(unique(df$abbrevs))) %>%
  ggplot(aes(x=lumMean, colour = subRename, y= after_stat(count)))+ 
  geom_density()+
  colScale+
  facet_wrap(~abbrevs, scales = 'free_y',ncol = 5)+
  theme_bw()+
#theme(legend.position = 'none')+
  scale_x_continuous(limits = c(0,0.7))+
xlab('Luminance')
dev.off()

png("output\\adaptive-landscape\\RvS\\facet-RG-density-weighted.png",
    type = 'cairo', units = 'mm', width = 215, height = 300, res = 300)
df %>%
  filter(as.integer(abbrevs) %in%
           1:length(unique(df$abbrevs))) %>%
  ggplot(aes(x=xCoord, colour = subRename, y= after_stat(count)))+ 
  geom_density()+
  colScale+
  facet_wrap(~abbrevs, scales = 'free_y',ncol = 5)+
  theme_bw()+
  #theme(legend.position = 'none')+
  #scale_x_continuous(limits = c(0,0.7))+
  xlab('R:G opponent')
dev.off()

png("output\\adaptive-landscape\\RvS\\facet-BY-density-weighted.png",
    type = 'cairo', units = 'mm', width = 215, height = 300, res = 300)
df %>%
  filter(as.integer(abbrevs) %in%
           1:length(unique(df$abbrevs))) %>%
  ggplot(aes(x=lumMean, colour = subRename, y= after_stat(count)))+ 
  geom_density()+
  colScale+
  facet_wrap(~abbrevs, scales = 'free_y',ncol = 5)+
  theme_bw()+
  #theme(legend.position = 'none')+
  #scale_x_continuous(limits = c(0,0.7))+
  xlab('B:Y opponent')
dev.off()
rm(x)

###** ridgeline plots of geology luminance ####
library(viridis)
library(ggridges)

png(file ='output//adaptive-landscape//geol-lum-ridgeplots.png',
    width=800, height=600, units = "px", res = 100)
ggplot(df, aes(x =  lumMean, y = geology, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1.) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_gradient(name = "Log luminance", low = "grey 30",high ="white") +
  theme_ridges(font_size = 13, grid = TRUE) + 
  xlab("Log luminance") + ylab("dominant geologies") +
  theme(legend.position = "none", axis.title.x = element_text(hjust = 0.1, face = "bold"),
        axis.title.y = element_text(hjust = 0.1, vjust = 3,face = "bold"))
dev.off()

# geologies
geol <- df[df$substrate == "b",]
geol <- geol[geol$geology != "graniteOrComplex"]

p <- ggplot(geol, aes(x = xCoord, y = yCoord, color = geology)) +
  geom_point(size = 1, stroke = 0.5, shape = 1) +  
  geom_point(x=0, y=0, shape = 3, color = "black", stroke = 1) +
  stat_ellipse(color = "black")  + theme_bw()

# TO DO: add marginal density plots
png(file= paste0('output//adaptive-landscape//geol-chrom.png'), 
    width=800, height=600, units = "px", res = 120)
p + facet_grid(rows = vars(geology)) +
  theme(axis.title=element_blank(),
        legend.position = "none") 
dev.off()


