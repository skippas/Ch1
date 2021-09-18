library(tidyverse)
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
# P and G are grouped with b!!!
df[df$substrate == "p" | df$substrate == "g", "substrate"] <- "b" 
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
df <- df %>% mutate(species = tolower(species))

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

rm(sampling,ls,dat,abbrevs, i, path,short,rowsub, excl, calcrete, dupROINames, 
   graniteOrComplex,quartz_quartzite, row_sub, sedimentary)

