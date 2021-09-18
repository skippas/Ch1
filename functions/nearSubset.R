# this function takes a grouped distance df as input, and returns a nearest 
# subset using either nearest percentages (qtiles) of observations
# or nearest n observations

# Need to be very careful with the group by term, can introduce errors that
# won't be noticed

nrst_subset <- function(datJnd, visfeature = "dL", qtiles = NULL, nNrst = NULL ){
  if(!is.null(qtiles)){
    jnds_ntiles <- datJnd %>% mutate(qtile = dplyr::ntile(!! sym(visfeature), qtiles))
    jnds_ntiles[jnds_ntiles$qtile %in% seq(1, 1, by = 1),]
  }else{
    datJnd %>% group_by(mspec.x) %>% slice_min(visfeature, n = nNrst)
  }
}

# Same as above function but with pre-filtering of dists
# TO DO: add return original DF functionality
nearSubset <- function(jnd_df = jnds_b, localOnly = T,acrossImage = F,
                       visfeature = "dL", qtiles = NULL, nNrst = NULL ){
    if(localOnly){
      jnd_df <- jnd_df[jnd_df$abbrevs.x == jnd_df$abbrevs.y,]
    }
    if(acrossImage == F){
      jnd_df <- jnd_df[jnd_df$mspec.x == jnd_df$mspec.y,]
      jnd_df <- jnd_df %>% group_by(abbrevs.y, abbrevs.x,mspec.x, mspec.y, roi.x) 
    }
    if(acrossImage == T){
      jnd_df <- jnd_df %>% group_by(abbrevs.y, abbrevs.x, mspec.x, mspec.y, roi.x) 
    }
    if(!is.null(nNrst) | !is.null(qtiles)){
      nrstSub <- nrst_subset(jnd_df, visfeature, qtiles, nNrst)
    }else{
      return(jnd_df)
    }
  }


