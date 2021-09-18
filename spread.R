source("loading-cleaning.R")
# spread 

# luminance MAD
# N.N. Taleb video: why everyone use SD? Weights outliers strongly.
source("functions//MAD.R")
lumSpread <- df %>% group_by(pop_spp, substrate) %>%
  summarise(MAD = MAD(lumMean))

boxplot(lumSpread$MAD~lumSpread$substrate)

# luminance MAD after nearest dists
source("nearest_qtile.R")
jnds_c <- colDistEffic(bg_subs = "c") 
jnds_c <- jnds_c[jnds_c$abbrevs.x == jnds_c$abbrevs.y,]
jnds_c <- jnds_c[jnds_c$mspec.x == jnds_c$mspec.y,]
# nearest 10% each substrate each image to each plant
jnds_grouped <- jnds_c %>% group_by(abbrevs.y, abbrevs.x,mspec.x, mspec.y, roi.x) 
nrstSub_c <- nrst_subset(datJnd = jnds_grouped, visfeature = "dL", nNrst = 1) 

jnds_b <- colDistEffic(bg_subs = "b")
jnds_b <- jnds_b[jnds_b$abbrevs.x == jnds_b$abbrevs.y,]
jnds_b <- jnds_b[jnds_b$mspec.x == jnds_b$mspec.y,]
jnds_grouped <- jnds_b %>% group_by(abbrevs.y, abbrevs.x,mspec.x, mspec.y, roi.x) 
nrstSub_b <- nrst_subset(datJnd = jnds_grouped, visfeature = "dL", nNrst = 1)

nrstSub_b$bg_substrate <- "b"
nrstSub_c$bg_substrate <- "c"
nrstSub <- rbind(nrstSub_b, nrstSub_c)
indexDF <- nrstSub[, c("abbrevs.x", "mspec.x", "roi.y")]
nrstDF <- merge(df[, c("swMean", "mwMean", "lwMean","lumMean",
                                "abbrevs", "mspec","roi", "substrate")], indexDF,
                by.y = c("abbrevs.x","mspec.x","roi.y"),
                by.x = c("abbrevs","mspec", "roi"))

lumSpread <- nrstDF %>% group_by(abbrevs, substrate) %>%
  summarise(MAD = MAD(lumMean))

boxplot(lumSpread$MAD~lumSpread$substrate)

# Colour MADs
colSpread <- df %>% group_by(pop_spp, substrate) %>%
  summarise(xMAD = MAD(xCoord), yMAD = MAD(yCoord))
boxplot(colSpread$xMAD~colSpread$substrate)
boxplot(colSpread$yMAD~colSpread$substrate)


