naureen <- df[df$abbrevs == "Rf_Nar_ST",]
boxplot(naureen$lumMean ~ naureen$substrate)

temp <- medsLum[medsLum$abbrevs.y == "Rf_Nar_ST",]
temp <- df[df$abbrevs == "bg_Hal.O_SR" & df$substrate == "a",]
temp$substrate <- "aa"
naureen <- rbind(temp, naureen)
boxplot(naureen$lumMean ~ naureen$substrate)

meyeri <- df[df$abbrevs == "Ls2_Mey_QR",]
boxplot(meyeri$lumMean ~ meyeri$substrate)

temp <- medsLum[medsLum$abbrevs.y == "Ls2_Mey_QR",]
temp <- df[df$abbrevs == "bg_Hal.O_SR" & df$substrate == "a",]
temp$substrate <- "aa"
naureen <- rbind(temp, naureen)
boxplot(naureen$lumMean ~ naureen$substrate)

temp <- medsCol[medsCol$abbrevs.y == "Wh_Hkr_JI",]

temp <- df[df$abbrevs == "Pt_Les.V_LK",]
boxplot(temp$lumMean ~ temp$substrate)
temp <- temp[temp$substrate== "c",]

temp <- jnds[jnds$abbrevs.y == "Pt_Les.V_LK",]
temp <- temp[temp$abbrevs.x == "Pt_Les.V_LK",]
temp <- temp[temp$group.y == "c"]

temp<- medsLum[medsLum$abbrevs.y == "Ls1_Mey_PQ",]
