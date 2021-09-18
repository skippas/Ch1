# this function can produce different nearest distance subsets for a given 
# grouped df, using either nearest percentages or nearest n observations
# returns either the nearest of a set of quantiles, or the nearest n observ

nrst_subset <- function(datJnd, visfeature = "dL", qtiles = NULL, nNrst = NULL ){
  if(!is.null(qtiles)){
    datJnd %>% mutate(qtile = ntile(!! sym(visfeature), qtiles))
    jnds_c_dec <- jnds_ntiles[jnds_ntiles$qtile %in% seq(1, 1, by = 1),]
  }else{
    datJnd %>% slice_min(visfeature, n = nNrst)
  }
    }




