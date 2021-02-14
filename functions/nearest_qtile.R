nearestQtile <- function(df = jnds, visfeature = "dL", ntiles = 4,nearest_ntile = 1){
  df <- df %>% group_by(abbrevs.y, abbrevs.x,roi.x) %>% mutate(qtile = ntile(!! sym(visfeature), ntiles))
  df <- df[df$qtile %in% seq(1, nearest_ntile, by = 1),]
  meds <- df %>% group_by(abbrevs.y, abbrevs.x, roi.x) %>% summarise(median = median(!! sym(visfeature)))
  meds <- meds %>% group_by(abbrevs.y, abbrevs.x) %>% summarise(median = mean(median))
  return(meds)
}