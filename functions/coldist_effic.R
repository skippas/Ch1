library(gridExtra)
library(data.table)
library(pavo)


# Homemade trichromatic colour distance function
colDistEffic <- function(bg_subs = c("c","b"), sub2 = "a"){
  
  # these are your human RELATIVE e values for each cone. 
  noise_e <- c(s = 0.5656854 ,m = 0.1414214 ,l =0.1000000) 
  achromatic_e <- 0.12 

  df <- as.data.table(df)
  
  denominator <- sum((combn(noise_e, m = 2, FUN = prod))^2)
  lith <- df[df$substrate == sub2, c("lwMean","mwMean","swMean",
                                    "lumMean", "roi","mspec", "abbrevs", "group")]
  bg <- df[df$substrate %in% c(bg_subs), c("lwMean","mwMean","swMean","lumMean",
                                           "roi","mspec", "abbrevs", "group")]
  bg<-as.data.table(bg)
  lith <- as.data.table(lith)
  jnds <- left_join(x= lith,y= bg, by = character()) #using this instead of dt merge may be reason for error
  
  jnds[,dlw := log(lwMean.x)- log(lwMean.y)]
  jnds[,dmw := log(mwMean.x)- log(mwMean.y)]
  jnds[,dsw := log(swMean.x)- log(swMean.y)]
  
  jnds[,numerator := noise_e["s"]^2*(dlw - dmw)^2 +
         noise_e["m"]^2*(dlw - dsw)^2 +
         noise_e["l"]^2*(dsw - dmw)^2]
  
  jnds[,dS := sqrt(numerator / denominator)]
  
  # achromatic contrasts
  jnds[, dL := abs((log(lumMean.x) - log(lumMean.y)) / achromatic_e)]
  
  jnds[, grep("Mean",colnames(jnds)):=NULL]
  jnds[, c("dlw","dmw","dsw","numerator"):=NULL]
  colnames(jnds)
  
  return(jnds)
  
}

# cross ref my output w pavo

#df <- df[df$pop_spp == grep("Vrede",unique(df$pop_spp),ignore.case = T, value = T),]
# to get random sample each level
#rndid <- with(df, ave(X, substrate, FUN=function(x) {sample.int(length(x))}))
#df <- df[rndid==5,]
#outpAndre <- colDistEffic()

#outpPavo <- df
#outpPavo <- outpPavo[, c("swMean", "mwMean", "lwMean", "lumMean")]
#colnames(outpPavo) <- c("s", "m", "l","lum")
#outpPavo <- coldist(outpPavo, n = c(1,16,32), weber = .1, weber.achro = .12, 
#                    weber.ref = "longest",qcatch = "Qi", achromatic = TRUE)

#rm(outpAndre) ; rm(outpPavo) ; rm(rndid) ; rm(df)
