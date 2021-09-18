source("loading-cleaning.r")
unique(df$pop_spp)
# Khoeries Fulleri KL should be incl in quartz by not currently in dataset 
# due to small sample size? 
# Kangnas marm LM excl for same reason as above? Also should go in quartz
quartz_quartzite <-  c("kenhardt_fulleri_FG", "SidiBarani_fulleri_EF", 
"SidiBarani_verruculosa_EF", "boegoesberg_hallii.ochraceae_SR", "pofadder_olivacea_XXOP",
"raaswater_bromfeldii.menelli_UT", "namies_olivacea_XXON", "keimoes1_bromfeldii.insularis_XV",
"Pofadder1_olivacea_YX", "namies_fulleri_JK", "geselskap_fulleri_MN", "kangnas1_marmorata_LM",
"vioolsdrif_dinteri.brevis_FGH", "kangnas2_marmorata_OP", "lekkersing1_meyeri_PQ", "bitterfontein_divergens_HIJ",
"lekkersing2_meyeri_QR", "alexanderbay_herrei_RS", "strykloof_divergens_VX", "SteenkampskraalRd_divergens_XY",
"Pella_Federici_BCD", "aggeneys_olivacea_CDE")
sedimentary <- c("Toekoms1_Verruculosa_GF","Toekoms2_Verruculosa_HG","Liebeberg_hookeri_XXHL",
                 "Reynekespan_Hookeri_PO","Haasriver1_Oztenia_IJK")
# There's certainly a prominent secondary dark geology in some (eg. bavi salic ON)
# of these calcrete populations
calcrete <- c("platfontein_leslei.venterii_LK", "Boshoff_Lesleii.Venterii_ML",
              "Copperton_hallii_DE", "Baviaanskrans_Salicola_ON", "pofadder_fulleri_GH",
              "Upington_bromfeldii_VU", "Blesbergmine1_Marmorata.Elisiae_DEF",
              "uitkyk_fulleri_TU")
# the rest are complex or granite. Granite pops im not sure about yet,
# look like might be some other geologies in them
graniteOrComplex <- c(calcrete, sedimentary, quartz_quartzite)
graniteOrComplex <- unique(df[!(df$pop_spp %in% graniteOrComplex), "pop_spp"])

df$geology[df$pop_spp %in% graniteOrComplex] <- 'graniteOrComplex'
df$geology[df$pop_spp %in% quartz_quartzite] <- 'quartz_quartzite'
df$geology[df$pop_spp %in% calcrete] <- 'calcrete'
df$geology[df$pop_spp %in% sedimentary] <- 'sedimentary'

# ridgeline plots of geology brightness
library(viridis)
library(ggridges)
ggplot(df, aes(x =  lumMean, y = geology, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1.) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_fill_viridis(name = "Lum ?units?", option = "C") +
  labs(title = 'Luminance of dominant geologies') +
  theme_ridges(font_size = 13, grid = TRUE) + theme(axis.title.y = element_blank())

# density plots
dfGeoOrder <- df[order(df$geology),]
dfGeoOrder$abbrevs <- factor(dfGeoOrder$abbrevs, levels = unique(dfGeoOrder$abbrevs))

png(file='output//adaptive-landscape//lum-dens_by-site_geoOrder.png', width=1400, height=1800, units = "px", res = 150)
par(mfrow=c(12,5), mar=c(.3,.3,.3,.3))
for(i in 1:length(levels(dfGeoOrder$abbrevs))){
  plot(1,xlim = c(0,1), ylim = c(0,18),
       xlab = "logLum", ylab = "density",
       xaxt = 'n', yaxt = 'n')
  title(levels(dfGeoOrder$abbrevs)[i], adj = 0.9, line = -1, cex.main = 0.8)
  title(unique(dfGeoOrder[dfGeoOrder$abbrevs == levels(dfGeoOrder$abbrevs)[i],"geology"]), 
        adj = 0.9, line = -2, cex.main = 0.8)
  lines(density(dfGeoOrder[dfGeoOrder$abbrevs == levels(dfGeoOrder$abbrevs)[i] & 
                             dfGeoOrder$substrate == "b",]$logLum), 
        col = "black", lty = 1, lwd = 1)
}
dev.off() 
