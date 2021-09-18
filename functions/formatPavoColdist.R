# convert from dataframe to pavo coldist format
formatPavoColdist <- function(df){
  rowNamez <- paste0(df$abbrevs,df$mspec,df$roi)
  df <- df[, c("swMean", "mwMean", "lwMean","lumMean")]
  row.names(df) <- rowNamez
  colnames(df) <- c("s", "m", "l", "lum") 
  return(df)
}