# TO DO
# marmorata ELISAE (not elis-I-ae)
# kangnas (which one?) marmorata has only 3 images. Lump?

library(tidyverse)
library(data.table)
path <- "data\\2021-04-06\\"
df <- read.csv(paste0(path,"Image Analysis Results Nikon D7100 CoastalOpt 105mm D65 to Human D65.csv"))
df <- separate(df, Label, into = c("pop_spp", "label"), sep = "/")
df$pop_spp <- basename(df$pop_spp)
dupROINames <- df[which(!duplicated(df$label)),]$label # remove dups
df <- df[which(!duplicated(df$label)),]
df <- separate(df, label, into = c("vis_image","uv_image", "substrate"), sep = "_")
df <- unite(df, "mspec", c("vis_image","uv_image"))
df$roi <- df$substrate
df$substrate <- substr(df$substrate, 0,1)
unique(df$substrate) # empty space substrates?
df <- df[df$substrate %in% c("a","b","c"),]

# short names
abbrevs <- c("Ag_Oli_CDE","Ab_Her_RS","At_Loc_BA","Bv_Sal_ON", 
             "Bf_Div_HIJ","Bb_Mar_NO","Bbm_Mar_EFG","Bbm1_Mar_DEF",
             "bg_Hal.O_SR","Bh_Les.V_ML","Br_Ozt_UV", "Cop_Hal_DE", 
             "Dd_Com.W_FE", "Gs_Ful_MN", "Hr1_Ozt_IJK","Hr2_Ozt_JKL",
             "Ht_Auc_AB", "Ht_Auc_XXAH", "Kn_Mar_LM01", "Kn1_Mar_LM02", 
             "Kn2_Mar_OP","Km1_Brm.I_XV", "Kh_Ful_KL","Ke_Ful_FG",
             "Kz_Div.A_GHI", "Lb_Hkr_XXHL", "Ls1_Mey_PQ","Ls2_Mey_QR",
             "Md_Hkr.H_RQ", "Mj1_Com_DC", "Mj2_Com_ED", "Nm_Ful_JK",
             "Nm_Oli_XXON", "Om_Hk_CD", "Pe_Fed_BCD","Pt_Les.L_KJ",
             "Pt_Les.V_LK", "Pf_Ful_GH","Pf_Oli_XXOP", "Pf1_Oli_YX",
             "Pf.C_Ful_ABC", "Rw_Brm.M_UT", "Ry_Hkr_PO","Rf_Nar_ST",
             "Sb_Ful_EF01", "Sb_Ver_EF02", "Skk_Div_XY", "Stryk_Div_VX",
             "Stryd_Hkr.H_QP","Tk1_Ver_GF","Tk2_Ver_HG", "Uk_Ful_TU",
             "Up_Brm_VU", "Up_Ful_TS", "Vb_Loc_CB","Vd_Dnt.B_FGH",
             "Vr_Hal_IH", "Wh_Hkr_JI", "Wtn_Les_NM")
short <- c('Aucamp_Hopet_AB','Aucamp_Hopet_XXAH', 'Bromf_Uping_VU',
           'Bromf.Insul_Keim1_XV','Bromf.Menel_Raasw_UT','Compt_Matj1_DC',
           'Compt_Matj2_ED','Compt.Web_Droed_FE','Dint.Brev_Viools_FGH',
           'div.ameth_kotzer_GHI','Div_Bitterf_HIJ','Div_Skampsk_XY',
           'Div_Stry_VX','Feder_Pella_BCD','Fuller_Gesels_MN','Fuller_Kenh_FG',
           'Fuller_khoer_KL','Fuller_Namies_JK','Fuller_Pof_GH','Fuller_SidiB_EF01',
           'Fuller_Pof.C_ABC','Fuller_Uitkyk_TU','Fuller_Uping_TS','Hallii_Copp_DE',
           'Hallii_Vrede_IH','Hallii.Ochr_Boeg_SR','Herrei_Alexb_RS',
           'Les.Vent_Plat_LK', 'Hook_Lieb_XXHL','Hook_Omdr_CD','Hook_Reyn_PO',
           'Hook_Whoe_JI', 'Hook.Hk_Mary_RQ','Hook.Hk_Stry_QP','Les.Vent_Boshf_ML',
           'Les.Les_Plat_KJ','Les_Wton_NM','Local_Attie_BA','Local_Verseb_CB',
           'Marm_Blesb_NO','Marm_Blesbm_EFG','Marm_Kang_LM01','Marm_Kang1_LM02',
           'Marm_Kang2_OP','Marm.El_Blesbm1_DEF','Meyeri_Lsing1_PQ','Meyeri_Lsing2_QR',
           'Naareen_Rooif_ST', 'Oliv_Aggen_CDE', 'Oliv_Namies_XXON','Oliv_Pof1_YX',
           'Oliv_Pof_XXOP', 'Ozten_Brakf_UV','Ozten_Haasr1_IJK',
           'Ozten_Haasr2_JKL','Salic_Bavi_ON','Verru_Sidib_EF02',
           'Verru_Toek1_GF','Verru_Toek2_HG')
df[df$pop_spp == "namies_olivacea_x","pop_spp"] <- "namies_olivacea_XXON"
df[df$pop_spp == "pofadder_olivacea_x","pop_spp"] <- "pofadder_olivacea_XXOP"
df[df$pop_spp == "Liebeberg_hookeri_x","pop_spp"] <- "Liebeberg_hookeri_XXHL"
df[df$pop_spp == "Hopetown_aucampiae_x","pop_spp"] <- "Hopetown_aucampiae_XXAH"

df$code = sapply(strsplit(df$pop_spp, "_"), function(x) x[3])
abbrevs <- as.data.frame(abbrevs)
abbrevs$code <- sapply(strsplit(abbrevs$abbrevs, "_"), function(x) x[3])

df <- merge(df, abbrevs, by = "code", all.x = TRUE)
unique(df$pop_spp)

# GEOLOGIES 

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

# take care of the multi-species sites
df[grep("SidiBarani_fulleri_EF", df$pop_spp),"abbrevs"] <- 
  abbrevs[grep("Sb_Ful_EF01", abbrevs$abbrevs),"abbrevs"]
df[grep("SidiBarani_verruculosa_EF", df$pop_spp),"abbrevs"] <- 
  abbrevs[grep("Sb_Ver_EF02", abbrevs$abbrevs),"abbrevs"]
df <- df[,!(colnames(df) %in% "code")]

df[grep("kangnas_marmorata_LM", df$pop_spp),"abbrevs"] <- 
  abbrevs[grep("Kn_Mar_LM01", abbrevs$abbrevs),"abbrevs"]
df[grep("kangnas1_marmorata_LM", df$pop_spp),"abbrevs"] <- 
  abbrevs[grep("Kn1_Mar_LM02", abbrevs$abbrevs),"abbrevs"]
df <- df[,!(colnames(df) %in% "code")]

df$group <- "bg"
df$group[df$substrate == "a"] <- "a"

which(is.na(df$abbrevs))
df <- na.omit(df)
df <- df[!rowSums(df < 0), ]
df <- df[df$substrate %in% c("a","b","c"),]
df$pop_spp <- as.factor(df$pop_spp)
df$abbrevs <- as.factor(df$abbrevs)

df$species <- stringr::str_match(df$pop_spp, "_\\s*(.*?)\\s*_")[,2]
df$population <- sapply(str_split(df$pop_spp, "_",  n = 2), `[`, 1)
df <- df %>% mutate(species = tolower(species), population = tolower(population))

# !!!!!!!!!!!! UNRESOLVED ISSUE !!!!!!!!!!!!!!
# kangnas 2 marmorate OP 12 looks like divergens

# Pops with low sample sizes
# Pops sampled twice, cf Allan
sampling <- df %>% group_by(pop_spp) %>% summarise(nImages = n_distinct(mspec))
# excl kangnas marmorata & khoeries fulleri -
# low sample size 3 & 4 images respectively
df <- df[!df$pop_spp %in% c("kangnas_marmorata_LM","khoeries_fulleri_KL"),]

# to drop unused levels 
df[] <- lapply(df, function(x) if(is.factor(x)) factor(x) else x)

# pesky 0's
row_sub <- apply(df[, c("swMean", "mwMean", "lwMean")], 
                1, function(row) all(row > 0.00001 ))
df <- df[row_sub,]

# Remove colour oultiers
source("functions//RNL_colspace_tri.R")
df <- RNL_colspace_tri(df = df)
df <- df %>% 
  mutate(lumCoord = log(lumMean)/0.05)
plot(df$xCoord, df$yCoord)
abline(h = c(1, -5.5), v =(8))
colOutliers <- df[xCoord > 8 | yCoord >1 | yCoord < -5.5,]
df <- df[xCoord < 8 & yCoord < 1 & yCoord > -5.5,]
dev.off()

# check with Allan
df[df$species == "lesleii.venterii", "species"] <- "leslei.venterii"
df[df$species == "oztenia", "species"] <- "otzeniana"

# Lump subspecies
df <- df %>% 
  mutate(species = fct_relevel(species, sort)) %>%
  mutate(collapseSpp =
           fct_collapse(species, 
                        bromfieldiiSpp = c("bromfeldii", "bromfeldii.insularis","bromfeldii.menelli"),
                        comptoniiSpp = c("comptonii", "comptonii.webrii"),
                        divergensSpp = c("divergens", "divergens.amethystina"),
                        halliiSpp = c("hallii", "hallii.ochraceae"),
                        hookeriSpp = c("hookeri", "hookeri.hookeri"),
                        leslieiSpp = c("leslei", "leslei.leslei", "leslei.venterii"),
                        marmorataSpp = c("marmorata", "marmorata.elisiae"),
                        dinteriSpp = c('federici','dinteri.brevis'))) %>%
  group_by(collapseSpp) %>% 
  mutate(collapseSpp = case_when(
    n_distinct(abbrevs) >1 ~ as.character(collapseSpp),
    n_distinct(abbrevs) <2 ~ "singlePopSpp")) %>%
  mutate(collapseSpp = as_factor(collapseSpp),
    collapseSpp = fct_recode(collapseSpp,
                                      bromf. = 'bromfieldiiSpp',
                                      compt. = 'comptoniiSpp',
                                      verruc. = 'verruculosa',
                                      dint. = 'dinteriSpp',
                                      singlePops = 'singlePopSpp'))

# put in a table with abbrev, pop, species, sampling of subs..
# TO DO: add in the mspec counts. Difficult because summarise causes to lose. 
summaryTbl <- df %>% arrange(species) %>%
  group_by(species, population,abbrevs, substrate) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = substrate, values_from = n) %>%
  rename("abbreviation" = "abbrevs", "Lithops ROI count" = "a",
         "Rock ROI count" = "b", "Soil ROI count" = "c")
write.csv(summaryTbl, quote = F)

# Read in Allan geology categories
alanGeol <- read.csv("C:/Users/User/OneDrive - Stellenbosch University/Masters/Dominant geologies/Geology-descriptions AGE.csv")
alanGeol <- alanGeol %>% unite(pop_spp, c("population", "species", "code"), sep = "_")
alanGeol[alanGeol$pop_spp == "Hopetown_aucampiae_XXAH", "allan.categories"] <- 
  alanGeol[alanGeol$pop_spp == "Hopetown_aucampiae_AB", "allan.categories"][1] # assume same geol (same site, diff year?)
df$pop_spp <- gsub("Oztenia", "otzeniana", df$pop_spp) # correct scientific name
df <- merge(df, alanGeol[, c("pop_spp","allan.categories")], all.x = T) 

df <- df %>%
  mutate(allan_geol = 
           fct_collapse(allan.categories,
                        sedimentary = c("sedimentary"),
                        paleMet = c("gneiss (pale meta)",
                                    "quartzite (pale meta)",
                                    "granite (pale meta)"), 
                        quartz = c("quartz"),
                        calcrete = c("calcrete"),
                        darkIgnMet = c("gneiss (dark igneous)",
                                       "dolerite (dark igneous)"),
                        mixed = c("mixed / dark igneous",
                                  "mixed / sedimentary", 
                                  "gneiss (pale meta) / mixed",
                                  "sedimentary / mixed")),
         geol_abbrevs = fct_recode(allan_geol,
                                   SED = 'sedimentary',
                                   PMET = 'paleMet',
                                   QTZ = 'quartz',
                                   CAL = 'calcrete',
                                   MIX = 'mixed',
                                   DIGMET = 'darkIgnMet')) %>%
  select(-c("allan.categories"))  
  
  

df[df$substrate == 'a', 'subRename'] <- 'Lithops'
df[df$substrate == 'b', 'subRename'] <- 'Rock'
df[df$substrate == 'c', 'subRename'] <- 'Soil'


#rm(list=setdiff(ls(), "df"))
rm(sampling,ls,dat,abbrevs, i, path,short,rowsub, excl,
   sedimentary, short, row_sub,
   path, quartz_quartzite, graniteOrComplex, calcrete,
   dupROINames,RNL_colspace_tri,
   alanGeol, colOutliers, summaryTbl)

for(i in 1:10){
  dev.off()
}

# double check JNDs / coordinates versus MICA and PAVO ####

#  euclidean <- function(a, b) sqrt(sum((a - b)^2))
# euclidean(c(1.6,-2.14), c(1.11,-1.78))
# euclidean(as.numeric(gMeanChrom[113,3:4]), as.numeric(gMeanChrom[2,3:4]))
# 
# 
# x<-df %>% filter(pop_spp == grep('Droedam', unique(df$pop_spp), value = T),
#                  mspec == '0317_0321', roi == 'a1' | roi == 'b1',
#                  substrate == 'a' | substrate == 'b')
# plot(x$xCoord, x$yCoord)
# x1 <-as.numeric(x[1, c('xCoord', 'yCoord')])
# x2 <- as.numeric(x[2, c('xCoord', 'yCoord')])
# euclidean(x1, x2)
# 
# checkgMeanDist<-df %>% filter(abbrevs == unique(df$abbrevs)[2]) %>%
#   group_by(abbrevs,mspec, substrate) %>%
#   summarise_at(c('swMean', 'mwMean', 'lwMean', 'lumMean'), gmean) %>%
#   group_by(abbrevs, substrate) %>%
#   summarise_at(c('swMean', 'mwMean', 'lwMean', 'lumMean'), gmean)
# 
# checkMeanCoord<-df %>% filter(abbrevs == unique(df$abbrevs)[2]) %>%
#   group_by(abbrevs,mspec, substrate) %>%
#   summarise_at(c('xCoord', 'yCoord'), mean) %>%
#   group_by(abbrevs, substrate) %>%
#   summarise_at(c('xCoord', 'yCoord'), mean)
# 
# checkgMeanDist <- checkgMeanDist %>%
#   select('abbrevs','substrate', 'swMean', 'mwMean','lwMean','lumMean') %>%
#   rename(s = 'swMean', m = 'mwMean', l = 'lwMean',  lum = 'lumMean') %>%
#   mutate(substrate = make.unique(substrate)) %>%
#   as.data.frame()
# 
# rownames(checkgMeanDist) <- checkgMeanDist$substrate
# checkgMeanDist$substrate <- NULL
# checkgMeanDist$abbrevs <- NULL
# cntrst <- substring(rownames(checkgMeanDist),1,1)
# 
# pavoDists <- bootcoldist(checkgMeanDist,
#                          by = cntrst,
#                          n = c(1, 16, 32),
#                          weber = 0.05,
#                          weber.achro = 0.1,
#                          achromatic = T,
#                          qcatch = 'Qi',
#                          weber.ref = 'longest'
# )
# 
# euclidean(as.numeric(checkMeanCoord[1,3:4]),
#           as.numeric(checkMeanCoord[2, 3:4]))
# 
# 
# library(pavo)
# coldist()
# 
# sqrt(32/49)*0.02
# 
# 0.01616244/ sqrt(16/49)

# Check coldist_effic function
# 
# efficValues <- nrstSub_bc[1:3,]
# efficValues <- efficValues %>% select(-c(visInfo, distance))
# semi_join(jnds_b, efficValues)
# 
# dat <- df %>%
#   filter(abbrevs == unique(df$abbrevs)[1]) 
# 
# dat <-dat %>% 
#   filter( (roi == 'a1' & mspec == '3859_3862') |
#            (roi == 'b3' & mspec == '3886_3888') |
#             (roi == 'b5' & mspec == '3880_3883')) 
# 
# dat <- dat %>% ungroup() %>%
#   select('abbrevs','substrate', 'swMean', 'mwMean','lwMean','lumMean') %>%
#   rename(s = 'swMean', m = 'mwMean', l = 'lwMean',  lum = 'lumMean') %>%
#   mutate(substrate = make.unique(substrate)) %>%
#   as.data.frame()
# 
# rownames(dat) <- dat$substrate
# dat$substrate <- NULL
# dat$abbrevs <- NULL
# cntrst <- substring(rownames(dat),1,1)
# 
# pavoDists <- coldist(dat,
#                          n = c(1, 16, 32),
#                          weber = 0.05,
#                          weber.achro = 0.05,
#                          achromatic = T,
#                          qcatch = 'Qi',
#                          weber.ref = 'longest'
# )
# coldis
# 
# 

library(ggpmisc)
ggplot(df, aes(mwMean, lwMean))+
  geom_point()+
  stat_poly_eq(
    formula = y~x,aes(label = paste(..eq.label..,..rr.label..,
                                    ..p.value.label.., sep = "~~~~")),
                             parse = TRUE)
