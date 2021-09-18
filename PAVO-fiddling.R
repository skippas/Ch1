# README 
# Double checking JND distances of scripts that I modified
# Plotting Lithops in interesting spaces
# Figuring out what the values are after visual modelling
## How do values I get from MICA compare to PAVO. 
## What would you do in PAVO to get same values?
## Pretty sure I dont exactly understand what is happening with the illuminant -
## and why you end up with the values you do, eg. 0.63. I think they are -
## proportions of the illuminant that the receptor catches.

source('loading-cleaning.R')
library(pavo)

set.seed(42)
lithops <- df %>% filter(substrate == "a") %>%
              group_by(abbrevs, mspec) %>% sample_n(size = 1) %>%
  group_by(abbrevs) %>%
  summarise(across(contains('Mean'),mean, names = 'mean_{.col}')) 

# Maxwell triangles
lithops <- as.data.frame(lithops[, c('abbrevs','swMean', 'mwMean', 'lwMean')])
rownames(lithops) <- lithops[,'abbrevs']
lithops[,1] <- NULL
names(lithops) <- c('s','m','l')

tri.lithops <- colspace(lithops, space = 'tri')
plot(tri.lithops, bg =)
plot(tri.lithops$x, tri.lithops$y)


# pavo datasets

specs <- readRDS(system.file("extdata/specsdata.rds",
                             package = "pavo"
))
mspecs <- aggspec(specs, by = 3, FUN = mean)
spp <- gsub("\\.[0-9].*$", "", names(mspecs))[-1]
sppspec <- aggspec(mspecs, by = spp, FUN = mean)
spec.sm <- procspec(sppspec, opt = "smooth", span = 0.2)

data(sppspec)
vismod1 <- vismodel(flowers,
                    visual = "avg.uv", achromatic = "bt.dc",
                    illum = "D65", relative = FALSE
)

data(flowers)
vis.flowers <- vismodel(flowers,
                        visual = "apis",
                        #illum = 'D65',
                        scale = 10000,
                        qcatch= 'fi',
                        #vonkries = T
)

tri.flowers <- colspace(vis.flowers, space = "tri")

plot(tri.flowers, pch = 21, bg = spec2rgb(flowers))


head(tri.flowers)

getAnywhere(vismodel)
