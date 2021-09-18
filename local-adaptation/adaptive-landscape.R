library(tidyverse)
source("loading-cleaning.R")

n <- df %>% group_by(abbrevs, substrate) %>% summarise(n = n())

theme_custom <- theme_bw()
theme_custom <- theme_update(panel.background = element_rect(fill = "white"),
                             panel.grid.major = element_line(size = .5),
                             panel.grid.minor = element_blank())
theme_set(theme_custom) 

# colour scatter of rock pops
rocks <- df %>% filter(substrate == "b")
sampsize <- rocks %>% group_by(abbrevs) %>% summarise(n = n())

ggplot(rocks, aes(xCoord, yCoord)) + geom_point() +
  facet_wrap(collapseSpp~abbrevs)

# boxplots each channel ####
p <- ggplot(data = df %>% group_by(abbrevs,substrate) %>% mutate(n = n()) %>% 
              mutate(label = paste0(substrate,'\nN = ',n)) %>%
              pivot_longer(cols = c(lwMean, mwMean, swMean),
                           names_to = "channel", values_to = "refl") %>%
              ungroup() %>% arrange(species) %>%
              mutate(abbrevs = factor(abbrevs, levels = unique(abbrevs))),
      aes(x = factor(label), y = refl, fill = channel)) + geom_boxplot()
p <- p +facet_wrap(vars(abbrevs),scales = "free_x")
png("output//adaptive-landscape//channel-reflectance-bySubstrate.png", type = 'cairo',
    width = 40, height = 60, units = "cm", res = 300, pointsize = 6) ; p; dev.off()

# boxplots of luminance estimates only rock and lithops ####
p <- ggplot(data = df %>% filter(substrate != "c"), aes(y = abbrevs, x = lumMean, fill = substrate)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_grid(collapseSpp~., scales = "free", space = "free") +
  theme_bw()
png("output//adaptive-landscape//luminance-bySubstrate.png", type = 'cairo',
    width = 200, height = 300, units = "mm", res = 300, pointsize = 6) ; p; dev.off()

# violin plots to assess bimod ####
p <- ggplot(data = df %>% group_by(abbrevs,substrate) %>% mutate(n = n()) %>% 
              mutate(label = paste0(substrate,'\nN = ',n)) %>%
              ungroup() %>% arrange(species) %>%
              mutate(abbrevs = factor(abbrevs, levels = unique(abbrevs))),
            aes(x = factor(label), y = lwMean)) +
  geom_violin(fill = "grey72", color = NA, scale = "count",bw = .05)
p <- p +facet_wrap(vars(abbrevs),scales = "free_x") + theme_minimal()
png("output//adaptive-landscape//lwReflectance-bySubstrate-violin.png", type = 'cairo',
    width = 40, height = 60, units = "cm", res = 300, pointsize = 6) ; p; dev.off()

library(ggdist)
p <- ggplot(data = df %>% group_by(abbrevs,substrate) %>% mutate(n = n()) %>% 
              mutate(label = paste0(substrate,'\nN = ',n)) %>%
              ungroup() %>% arrange(species) %>%
              mutate(abbrevs = factor(abbrevs, levels = unique(abbrevs))),
            aes(x = factor(label), y = lwMean)) +
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  coord_cartesian(xlim = c(1.2, NA), clip = "off")
p <- p +facet_wrap(vars(abbrevs),scales = "free_x") + theme_minimal()
png("output//adaptive-landscape//lwReflectance-bySubstrate-raincloud.png", type = 'cairo',
    width = 40, height = 60, units = "cm", res = 300, pointsize = 6) ; p; dev.off()

# Lithops brightness and colour raw data by species with stat tests 
ggplot(df %>% filter(substrate == "a"), aes(x = abbrevs, y = lumMean)) +
  geom_boxplot() +
  facet_grid(.~collapseSpp, scales = "free", space = "free") +   # funny, the notation in facet grid NB, vars() gives weird result
  scale_x_discrete(guide = guide_axis(angle = 45))

# add dummy var to reset colors in each facet
pdata <- df %>% group_by(collapseSpp) %>%
  mutate(colr = as.integer(factor(abbrevs)))
ggplot(pdata %>% filter(substrate == "a"), aes(x = xCoord, y = yCoord, color = as.factor(colr))) +
  geom_point() + facet_wrap(.~collapseSpp, scales = "free") +
  theme(legend.position = "none")
  
# Stats
manovaResLum <- df %>% ungroup() %>% filter(substrate == "a") %>% 
  nest_by(collapseSpp) %>% 
  mutate(Model = list(manova(cbind(xCoord, yCoord) ~ abbrevs, data = data))) %>%
  rowwise() %>% 
  mutate(ModelSummary = list(broom::tidy(Model, "Pillai"))) %>%
  select(-c(Model)) %>% unnest(ModelSummary)

anovaResCol <- df %>% ungroup() %>% filter(substrate == "a") %>% 
  nest_by(collapseSpp) %>% 
  mutate(Model = list(aov(lumMean ~ abbrevs, data = data))) %>% 
  rowwise() %>%
  mutate(ModelSummary = list(broom::tidy(Model))) %>%
  select(-c(Model)) %>% unnest(ModelSummary)

write.csv(anovaRes, quote = F)
write.csv(manovaRes, quote = F)

# Geology visual properties ####

#* SumStats ####
rocks <- df %>% filter(substrate == 'b')

rocks %>% 
  group_by(allan_geol) %>%
  summarise(sdxCoord = sd(xCoord),
            sdyCoord = sd(yCoord),
            sdLum = sd(lumMean)) %>%
  transmute(allan_geol,
            sdChroma = rowMeans(select(., c(sdyCoord, sdxCoord))))
  

#* Ridgeline plots ####
library(ggridges)

nLabels <- df %>% 
  filter(substrate == 'b') %>%
  group_by(allan_geol) %>%
  summarise(nPops = n_distinct(abbrevs),
            nROIs = n()) %>%
  mutate(xPos = 0.75,
         label = paste('nPops = ', nPops,'  ', 'nROIs = ', nROIs))

geolLabels <- c(paleMet = 'pale metamorphic', quartz = 'quartz',sedimentary = 'sedimentary',
  calcrete = 'calcrete', darkIgnMet = 'dark igneous \nor metamorphic', 
  mixed = 'mixed')

lumRidgePlot <- df %>% 
  filter(substrate == 'b') %>% 
  mutate(allan_geol =
           fct_reorder(allan_geol, lumMean, median, .desc = T)) %>% 
ggplot(aes(x = lumMean, y = allan_geol))+
  stat_density_ridges(aes(rel_min_height = 0.01),
                      quantile_lines = TRUE, quantiles = 2)+
  geom_text(data = nLabels,
            aes(label = label, x=xPos), vjust = -.8, size = 3)+
  ylab('')+
  xlab('Luminance value')+
  scale_y_discrete(labels = geolLabels)+
  scale_x_continuous(limits = c(0,.95))+
  theme_ridges()

RGRidgePlot <- df %>% 
  filter(substrate == 'b') %>% 
  mutate(allan_geol =
           fct_reorder(allan_geol, xCoord, median, .desc = T)) %>% 
  ggplot(aes(x = xCoord, y = allan_geol))+
  stat_density_ridges(aes(rel_min_height = 0.015),
                      quantile_lines = TRUE, quantiles = 2)+
  ylab('Lithology categories')+
  xlab('R:G opponent value')+
  geolLabels+
  scale_x_continuous(limits = c(0,2.5))+
  theme_ridges()

BYRidgePlot <- df %>% 
  filter(substrate == 'b') %>% 
  mutate(allan_geol =
           fct_reorder(allan_geol, yCoord, median, .desc = T)) %>% 
  ggplot(aes(x = yCoord, y = allan_geol))+
  stat_density_ridges(aes(rel_min_height = 0.015),
                      quantile_lines = TRUE, quantiles = 2)+
  ylab('')+
  xlab('B:Y opponent value')+
  geolLabels+
  scale_x_continuous(limits = c(-2,.3))+
  theme_ridges()

library(patchwork)

png("output//adaptive-landscape//geologies-ridges-composite.png", type = 'cairo',
    width = 215, height = 300, units = "mm", res = 300, pointsize = 6)
lumRidgePlot + RGRidgePlot + BYRidgePlot + plot_layout(ncol = 1)
dev.off()

#* geology 2D chroma densities / chroma scatterplots ####
p <- ggplot(rocks, aes(xCoord, yCoord)) + geom_point(alpha = .2) + 
  facet_wrap(allan_geol~abbrevs, ncol = 8)
png("output//adaptive-landscape//2d-colour-scatter.png",type = 'cairo',
    width = 200, height = 150, units = "mm", res = 300, pointsize = 6)
p
dev.off()

twoDContP <- ggplot(rocks, aes(xCoord, yCoord))+
  geom_density_2d_filled()+
  facet_wrap(vars(allan_geol),
             labeller = labeller(allan_geol = geolLabels))+ 
  xlim(0,2)+
  ylim(-2,0)+
  xlab('R:G opponent value')+
  ylab('B:Y opponent value')
png("output//adaptive-landscape//geologies-2d-contour-plot.png",type = 'cairo',
    width = 200, height = 150, units = "mm", res = 300, pointsize = 6)
twoDContP
dev.off()

library("ggpointdensity")
p <- ggplot(rocks, aes(xCoord, yCoord)) +
  geom_pointdensity() +
  facet_wrap(vars(allan_geol)) + scale_color_viridis_c(option = "plasma")
png("output//adaptive-landscape//geologies-pointdensity-plot.png",type = 'cairo',
    width = 200, height = 150, units = "mm", res = 300, pointsize = 6); p; dev.off()

ggplot(rocks %>% pivot_longer(names_to = "oppMech", values_to = "stim",
                               cols = c("xCoord", "yCoord")),
       aes(x = stim)) + geom_density() +
  facet_grid(allan_geol~oppMech)

p <- ggplot(rocks, aes(xCoord, yCoord)) +
  geom_pointdensity() +
  facet_wrap(collapseSpp ~ abbrevs) +scale_color_viridis_c(option = "plasma")
png("output//adaptive-landscape//pops-pointdensity-plot.png",type = 'cairo',
    width = 200, height = 300, units = "mm", res = 300, pointsize = 6); p; dev.off()

rocks <- rocks %>% group_by(collapseSpp, allan_geol) %>%
  summarise(popReps = n_distinct(abbrevs))
p <- ggplot(rocks, aes(x = allan_geol, y = collapseSpp,
              size = popReps, colour = popReps)) + geom_point() +
  scale_color_continuous(guide = "legend")  + theme_bw()
png("output//adaptive-landscape//species-geology-associations.png",type = 'cairo',
    width = 200, height = 200, units = "mm", res = 300, pointsize = 6); p; dev.off()


