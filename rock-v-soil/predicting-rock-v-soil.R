# What predicts which substrate is matched?

# calculate the necessary distances, then ask:

# 1. The 'surface area' of that substrate?
## after normalising by the distance betw the substrates?
## and after removing pops which match neither substrate
## after using only pops which have large dist betw rock and soil (minimizing noise)

# The brightness of rocks (ie. brighter not matched)
# Evol constraints ie. species? 
# Spread?

# repeat for colour and luminance

# TO DO: balance observations between images (like in nearest subsets) in gmean dists betw all obs

# SETUP ####
detach("package:urbnthemes", unload=TRUE)
detach("package:tidyverse", unload=TRUE)
library(tidyverse)
geom_aes_defaults()
update_geom_defaults("point", list(colour = "black"))

source("loading-cleaning.R")
source("functions//gMeanDists.R")
source("functions//nearSubset.R")
source("functions//distToOrigDF.R")
source("functions//coldist_effic.R")

# gmean dists between all obs ####
gMeanDistsAll <- df %>% 
  group_by(abbrevs,substrate, mspec) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), mean) %>%
  gMeanDists(df = ., combine_bg = F, substrates = c("b","c")) %>% 
  filter(comparison == "local") %>%
  pivot_wider(values_from = c("dL","dS"),
              names_from = c("patch1", "patch2"))

# gmeanDists betw nearest subsets ####
### 1 lithops per image not reproducible. see https://stackoverflow.com/questions/23831711/sample-n-random-rows-per-group-in-a-dataframe
set.seed(42)
oneLith <- df %>% filter(substrate == "a") %>% group_by(abbrevs, mspec) %>% 
  sample_n(size = 1) 
df <- rbind(df[df$substrate != "a",], oneLith)

# nearest 10% each substrate 
jnds_c <- colDistEffic(bg_subs = "c") 
jnds_b <- colDistEffic(bg_subs = "b")
nrstSub_c <- nearSubset(jnds_c, localOnly = T, acrossImage = T, 
                        visfeature = "dL", qtiles =  10) 
nrstSub_b <- nearSubset(jnds_b, localOnly = T, acrossImage = T, 
                        visfeature = "dL", qtiles =  10) 
nrstSub <- rbind(nrstSub_b, nrstSub_c)

### merge: This is a long way round. Will probably give similar answer to 
### just using the mean of distances across images. And that would be more consistent.
nrDistDF <- distToOrigDF(distDF = nrstSub, origDF = df)
nrDistDF <- nrDistDF %>% 
  group_by(abbrevs,substrate, mspec) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), mean) 
### calculate gmean dists
gMeanDistsNear <- gMeanDists(df = nrDistDF,
                             combine_bg = F,substrates = c("b","c")) %>%
  filter(comparison == "local") %>%
  pivot_wider(values_from = c("dL","dS"), names_from = c("patch1", "patch2"))

# Predicting distance from substrate frequency ####
###* Plot data #### 
source("rock-v-soil//rock-estimate.R")

props <- df %>% count(abbrevs, substrate)
props <- props %>% pivot_wider(names_from = substrate, values_from = n)
props$rock_cover <- (props$b / (props$b + props$c)*100)
props$soil_cover <- (props$c / (props$b + props$c)*100)
props <- arrange(props, as.factor(abbrevs))

propsAndDists <- merge(props, gMeanDistsAll, by.x = "abbrevs", by.y = "bgpop")
source('loading-cleaning.R')

###* Plots ####

###** Relship match and abundance scatterplot ####
library(ggpmisc)
plotF <- function(xvar = "soil_cover", yvar = "visfeature", xlab, ylab){
  ggplot(propsAndDists, aes(x = !! sym(xvar), y = !! sym(yvar))) +
    geom_point()+ 
    ylab(ylab)+
    xlab(xlab)+
    theme_classic()+
    geom_smooth(method = "lm", color = "black",
                formula = y~x, se =F, lwd = .5)+
    stat_poly_eq(formula = y~x,
                 aes(label = paste(..eq.label..,
                                   ..rr.label..,
                                   ..p.value.label.., sep = "~~~~")),
                 parse = TRUE)
}

yvar <- c("dL_a_c", "dL_a_b", "dS_a_c", "dS_a_b")
xvar <- paste0(c("soil", "rock", "soil", "rock"), "_cover")
xlabs <- paste(c("soil", "rock", "soil", "rock"), "surface area")
ylabs <- paste0("L",c("S lum", "R lum", "S chroma", "R chroma"),
                " dist")
plotList <- Map(plotF,xvar, yvar, xlabs, ylabs)

propsAndDists <- merge(props, gMeanDistsNear, by.x = "abbrevs", by.y = "bgpop")
plotListNrst <- Map(plotF, xvar, yvar, xlabs, ylabs)
png("output//predicting-rs//rs-meanDist-asfn-SurfA.png",
    units = 'mm', width = 215, height = 140, res = 300, type = "cairo")
grid.arrange(plotListNrst[[1]],plotList[[3]], plotListNrst[[2]], plotList[[4]])
dev.off()

###** RS dist also explains L- dist ####
gMeanDistsAll <- df %>% 
  group_by(abbrevs,substrate, mspec) %>% 
  summarise_at(vars(swMean, mwMean, lwMean, lumMean), mean) %>%
  gMeanDists(df = ., combine_bg = F, substrates = c("b","c")) %>% 
  filter(comparison == "local") %>%
  pivot_wider(values_from = c("dL","dS"),
              names_from = c("patch1", "patch2"))

propsAndDists <- merge(props, gMeanDistsAll, by.x = "abbrevs", by.y = "bgpop")

plotF <- function(xvar = "soil_cover", yvar = "visfeature",
                  xlab, ylab, yint){
  ggplot(propsAndDists, aes(x = !! sym(xvar), y = !! sym(yvar)))+
  geom_point()+
  geom_hline(yintercept = yint,lty =2)+
  geom_abline(lty = 2)+
  stat_poly_eq(formula = y~x,
               aes(label = paste(..eq.label.., ..rr.label..,
                                 ..p.value.label.., sep = "~~~")), 
               parse = TRUE)+
  theme_classic()+
  geom_smooth(method='lm', color = 'black', se = F,lwd = .5)+
  xlab(xlab)+
  ylab(ylab)
}

yvar <- c("dL_a_c", "dL_a_b", "dS_a_c", "dS_a_b")
xvar <- c("dL_b_c", "dL_b_c", "dS_b_c", "dS_b_c")
ylabs <- paste0('Mean (all points) L',c('S','R','S','R'),
                c(" lum", " lum", " chrom", " chrom"),
               ' dist')
xlabs <- paste('Mean (all points)',c("lum", "lum", "chrom", "chrom"),
               'RS dist')
yint <- c(1,1,0.1,0.1)
plotList <- Map(plotF,xvar, yvar, xlabs, ylabs, yint)

p <- plotList[[1]]
p2 <- plotList[[2]]
p3 <- plotList[[3]]
p4 <- plotList[[4]]
png("output\\predicting-rs\\RS-dist-pred-Lith-sub-dist-scatterplot.png",
    width = 215, height = 149, res = 300, type = "cairo", units = 'mm')
cowplot::plot_grid(p, p3,p2, p4, align = 'hv')
dev.off()

###** Specialisation metric ####
## Standardizing the distance between lithops and a substrate by the 
## distance between that substrate and the other backg substrate ('specialisation')

plotF(xvar = "rock_cover", yvar = "dL_a_c / dL_b_c" , xlab = "", ylab = "")
ggplot(propDist, aes(x = rock_cover, y = dL_a_b / dL_b_c)) +
  geom_point() + theme_bw() + geom_smooth(method = "lm") + 
  ylab("achro dists to soil") + xlab("% of surface area = soil")
dev.off()

sigBgDists <- propsAndDists[propsAndDists$dL_b_c > 5,]
ggplot(sigBgDists, aes(x = rock_cover, y = dL_a_c / dL_b_c)) +
  geom_point() + theme_bw() + geom_smooth(method = "lm") + 
  ylab("") + xlab("")
ggplot(sigBgDists, aes(x = soil_cover, y = dL_a_b / dL_b_c)) +
  geom_point() + theme_bw() + geom_smooth(method = "lm") + 
  ylab("") + xlab("")

propsAndDists$specialization <- propsAndDists$dL_a_c / propsAndDists$dL_b_c


library(urbnthemes)
theme_set(theme_minimal())
urbnthemes::set_urbn_defaults(style = "print")
###** background proportions stack plot ####
df<- df %>% filter(substrate != "a") %>% # this df is from rock estimate file
  mutate(abbrevs = fct_reorder(.f = abbrevs, 
                                .x = substrate,
                                .fun = function(.x) mean(.x == "b"),
                                .desc = TRUE),
         substrate = as_factor(substrate),
         substrate = relevel(substrate, "b"),
         subRename = fct_recode(substrate, rock = "b", soil = "c"))

png("output//adaptive-landscape//RvS//RS-prop-stackBar.png",
    type = 'cairo',width = 150, height = 180,
    units = "mm", res = 300, pointsize = 6)
ggplot(df)+
  geom_bar(mapping = aes(x = abbrevs, fill = subRename),
           position = position_fill(reverse = T))+
    geom_hline(aes(yintercept = .5), lty = "dotted")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     labels = scales::percent)+
  labs(x = "populations", y = NULL)+
    remove_ticks()+ guides(color = FALSE)+ coord_flip() 
dev.off()

###** The mean of background proportions  across pops ####
df %>%
  count(abbrevs, substrate) %>%
  pivot_wider(names_from = substrate, values_from = n) %>%
  mutate(propB = b / (b+c),
         propC = c / (b+c)) %>%
  summarise(mean(propB),
            mean(propC))

