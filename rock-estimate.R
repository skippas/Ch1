library(tidyverse)
path <- "data\\"
df <- read.csv(paste0(path,"Image Analysis Results Nikon D7100 CoastalOpt 105mm D65 to Human D65.csv"))
df <- separate(df, Label, into = c("pop_spp", "label"), sep = "/")
rowsub <- grep("whole",df$label)
df$label <- make.unique(df$label)
df <- separate(df, label, into = c("vis_image","uv_image", "substrate"), sep = "_")
df <- unite(df, "mspec", c("vis_image","uv_image"))
df$roi <- df$substrate
df$substrate <- substr(df$substrate, 0,1)
df[df$substrate %in% c("g","p"),"substrate"] <- "b" ## unique line!
df <- df[df$substrate %in% c("a", "b", "c"),] 
unique(df$pop_spp)
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
df$group[df$substrate == "a"] <- "lithops"

which(is.na(df$abbrevs))

# Pops with low sample sizes
# Pops sampled twice, cf Allan
sampling <- df %>% group_by(pop_spp) %>% summarise(nImages = n_distinct(mspec))
# excl kangnas marmorata & khoeries fulleri -
# low sample size 3 & 4 images respectively
df <- df[!df$pop_spp %in% c("kangnas_marmorata_LM","khoeries_fulleri_KL"),]

unique(df$abbrevs)
df <- na.omit(df)

# !!!!!!!!!!!! UNRESOLVED ISSUE !!!!!!!!!!!!!!
# kangnas 2 marmorata OP 12 looks like divergens

#rm(sampling,ls,dat,abbrevs, i, path,short,rowsub, excl)
