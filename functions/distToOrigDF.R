# Get the subset of observations in the original df which are present in 
# the nearest distance df
# careful, this probably only works for specific intra pop / image case

distToOrigDF <- function(distDF = nrstSub, origDF = df){
  xSubset <- merge(unique(origDF[, c("swMean", "mwMean", "lwMean","lumMean",
               "abbrevs", "mspec","roi", "substrate")]), 
        unique(distDF[, c("abbrevs.x", "mspec.x", "roi.x")]),
        by.y = c("abbrevs.x","mspec.x","roi.x"),
        by.x = c("abbrevs","mspec", "roi"))
  ySubset <- merge(unique(origDF[, c("swMean", "mwMean", "lwMean","lumMean",
                   "abbrevs", "mspec","roi", "substrate")]), 
        unique(distDF[, c("abbrevs.x", "mspec.x", "roi.y")]),
        by.y = c("abbrevs.x","mspec.x","roi.y"),
        by.x = c("abbrevs","mspec", "roi"))
  origDFSubset <- rbind(xSubset, ySubset)
  return(origDFSubset)
} 
