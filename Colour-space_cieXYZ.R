source("loading-cleaning.R")
library(pavo)

data("flowers")
vis.flowers <- vismodel(flowers, visual = 'cie10', illum = 'D65', vonkries = TRUE, relative = FALSE, achromatic = 'none')

colnames(df)[c(7,9,11)] <- c("X","Y","Z")
ls <- split(df, df$pop_spp)
ls <- lapply(ls, function(x) split(x, x$substrate))
ls <- lapply(ls, lapply, function(x) x %>% 
               select(X, Y, Z))
ls <- lapply(ls, lapply, colspace, space = "ciexyz", "cie10","d65")

pdf("output//cieXYZ-space.pdf")
par(mfrow = c(1,3))
for(i in 1:length(ls)){
  popspp <- names(ls[i])
  plot(ls[[i]][[1]], main = paste(popspp, names(ls[[i]])[1])) 
  plot(ls[[i]][[2]], main = paste(popspp, names(ls[[i]])[2])) 
  plot(ls[[i]][[3]], main = paste(popspp, names(ls[[i]])[3])) 
  
}
dev.off()

lithops <- df[df$substrate == "a" & df$pop_spp == "Hopetown_aucampiae_AB",]
colnames(lithops)[c(7,9,11)] <- c("X","Y","Z")
lithops <- lithops %>% select("X","Y","Z")
cieXYZ.lithops <- colspace(lithops, "ciexyz", "cie10")
plot(cieXYZ.lithops)
