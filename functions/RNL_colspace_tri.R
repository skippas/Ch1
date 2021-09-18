# Returns dataframe with colour space coordinates of ROIs
# Not applied to luminance yet (although might be trivial)
# see Renoult colour spaces in ecology for formula
RNL_colspace_tri <- function(df = df){
  v <- 0.05* sqrt(32/49) # noise in a single receptor
  wf <- c(s = v/sqrt(1/49) ,m = v/sqrt(16/49) ,l = v/sqrt(32/49)) # these are your human e values for each cone
  df <- as.data.table(df)
  
  # x coord (R:G Opp mech)
  df[,lwLog := log(lwMean)]
  df[,mwLog := log(mwMean)]
  df[,swLog := log(swMean)]
  
  opp1NoiseTerm <- sqrt(1/(wf["m"]^2 + wf["l"]^2))
  df[,xCoord := (lwLog - mwLog)*opp1NoiseTerm]
  
  # y coord 
  wfC <- combn(wf,2)
  opp2NoiseTerm <- unname(sqrt((wf["m"]^2 + wf["l"]^2) / 
                                 (prod(wf["s"],wf["m"])^2 +
                                    prod(wf["s"],wf["l"])^2 + 
                                    prod(wf["m"],wf["l"])^2)))
  df[,yCoord := (df[,lwLog*(wf['m']^2 / (wf['m']^2 + wf['l']^2))] + 
                   df[,mwLog*(wf['l']^2 / (wf['m']^2 + wf['l']^2))])]
  #df[,yCoord := (swLog*opp2NoiseTerm) - (yCoord*opp2NoiseTerm)]
  df[,yCoord := (swLog - yCoord)*opp2NoiseTerm]
  
  #df[,yCoord := yCoord*opp2NoiseTerm] I changed this term to above two, but it
  # wasn't actually the problem. Could change back.
  
  return(df)
}


