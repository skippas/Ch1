# Use this function to return a dataframe with the geometric mean distances
# between substrates
# This avoids problems with using the means of all distances
# Instead, the log catch of the receptors is found
# Note: This value is then exponentiated to reverse the log, because the log
## is later taken again in the coldist function.

# To do ####
# BUG: 
# when input dataframe has a single sub, but substrates argument specifies both
# get a different answer to when a single substrate is specified.
# if input df has both substrates and both substrates argument is specified, get
# the same answer to a substrate as an input df with one substrate and one substrate
# argument is specified


library(data.table)
library(progressr)
library(future.apply)
library(pavo)
library(gtools)

gMeanDists <- function(df = df, substrates = c("b", "c"), combine_bg = F){

# calculate gmeans
source("functions//gmean.R")

  dat <- df[df$substrate %in% c("a", substrates),]
  dat %>% group_by(abbrevs, substrate) -> dat
  if(combine_bg){
    dat %>% group_by(abbrevs, group) -> dat
  }
  dat %>% summarise_at(vars(swMean, mwMean, lwMean, lumMean), gmean) -> dat
  dat <- as.data.frame(dat)
  
# dists betw gmeans
  lBgCombs <- permutations(n = length(unique(as.character(df$abbrevs))),
                           r=2,
                           v = unique(as.character(df$abbrevs)),
                           repeats.allowed = T)
  lstCombs <- vector(mode = "list", length = length(lBgCombs[,1])) 
  names(lstCombs) <- paste(lBgCombs[,1], lBgCombs[,2], sep = "--")
  
  for(i in 1:length(lBgCombs[,1])){
    if(combine_bg){
      lstCombs[[i]] <- rbind(dat[dat$abbrevs == lBgCombs[i,1] & dat$group == "a",],
                             dat[dat$group == "bg" & dat$abbrevs == lBgCombs[i,2],])
    }
    else{
      lstCombs[[i]] <- rbind(dat[dat$abbrevs == lBgCombs[i,1] & dat$substrate == "a",],
                             dat[dat$substrate %in% substrates & dat$abbrevs == lBgCombs[i,2],])
    }
  }
  
  lstCombs <- lapply(lstCombs, function(x) {
    rownames(x) <- make.unique(x$substrate) ; x
  })
  lstCombs <- lapply(lstCombs, function(x) {
    x <- x[,-which(names(x) %in% c("abbrevs", "substrate"))] ; x
  })
  lstCombs <- lapply(lstCombs, function(x) {
    colnames(x) <- c("s", "m", "l","lum") ; x
  })
  
  lstDists <- lapply(lstCombs, coldist,
                     n = c(1,16,32), weber = .05, weber.achro = .1, 
                     weber.ref = "longest", qcatch = "Qi", achromatic = TRUE)
  
  dfDist <- rbindlist(lstDists, idcol= "names")
  dfDist <- separate(dfDist, names, into = c("lithpop", "bgpop"), sep = "--")
  dfDist[dfDist$lithpop == dfDist$bgpop,"comparison"] <- "local"
  dfDist[dfDist$lithpop != dfDist$bgpop,"comparison"] <- "foreign"
  
  return(dfDist)
}


# Double check JNDs versus PAVO's bootcoldist

# checkPavoVsGmeandist <-df %>% ungroup() %>%
#   filter(abbrevs == unique(df$abbrevs)[8]) %>%
#   select('abbrevs','substrate', 'swMean', 'mwMean','lwMean','lumMean') %>%
#   rename(s = 'swMean', m = 'mwMean', l = 'lwMean',  lum = 'lumMean') %>%
#   mutate(substrate = make.unique(substrate)) %>%
#   as.data.frame()
# rownames(checkPavoVsGmeandist) <- checkPavoVsGmeandist$substrate
# checkPavoVsGmeandist$substrate <- NULL
# checkPavoVsGmeandist$abbrevs <- NULL
# cntrst <- substring(rownames(checkPavoVsGmeandist),1,1)
# 
# pavoDists <- bootcoldist(checkPavoVsGmeandist,
#                          by = cntrst,
#                          n = c(1, 16, 32),
#                          weber = 0.1,
#                          weber.achro = 0.1,
#                          achromatic = T,
#                          qcatch = 'Qi', 
#                          weber.ref = 'longest'
# )

