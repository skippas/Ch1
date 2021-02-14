library(gridExtra)
library(data.table)
library(pavo)
source("Loading-cleaning.R")

colDistEffic <- function(bg_subs = c("c","b"), sub2 = "a"){
  
  # the following objects need to be renamed. 
  # calculations are to determine noise, too lazy to rename downstream.
  weberFrac <- c(.02,.02,.08) 
  relDens <- c(.33,.33,.33)
  snr_v <- weberFrac * sqrt(relDens)
  noise_e <- snr_v / sqrt(relDens)
  
  denominator <- sum((combn(noise_e, m = 2, FUN = prod))^2)
  
  lith <- df[df$substrate == sub2, c("lwMean","mwMean","swMean",
                                    "lumMean", "roi","mspec", "abbrevs", "group")]
  bg <- df[df$substrate %in% c(bg_subs), c("lwMean","mwMean","swMean","lumMean",
                                           "roi","mspec", "abbrevs", "group")]
  bg<-as.data.table(bg)
  lith <- as.data.table(lith)
  jnds <- left_join(x=lith,y=bg, by = character())
  
  jnds[,dlw := log(lwMean.x)- log(lwMean.y)]
  jnds[,dmw := log(mwMean.x)- log(mwMean.y)]
  jnds[,dsw := log(swMean.x)- log(swMean.y)]
  
  jnds[,numerator := noise_e[3]^2*(dlw - dmw)^2 +
         noise_e[2]^2*(dlw - dsw)^2 +
         noise_e[1]^2*(dsw - dmw)^2]
  
  jnds[,dS := sqrt(numerator / denominator)]
  
  # achromatic contrasts
  # Wi = channel noise (Vi) / sqrt(relative receptor density (Ni))
  # achr contrast (L) = receptor contrast (dF) / W
  jnds[, dL := abs((log(lumMean.x) - log(lumMean.y)) / .02)]
  
  jnds[, grep("Mean",colnames(jnds)):=NULL]
  jnds[, c("dlw","dmw","dsw","numerator"):=NULL]
  colnames(jnds)
  
  return(jnds)
  
}


