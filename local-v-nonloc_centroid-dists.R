df[df$substrate != "a", "sub"] <- "bg"
df[df$substrate == "a", "sub"] <- "l"

meds <- df %>% group_by(pop_spp, sub) %>% summarise(median = median(lumMean))
meds <- pivot_wider(meds, names_from = sub, values_from = median)
l <- meds[,c(1,3)]
bg <- meds[,c(1,2)]
jnds <- left_join(x=l,y=bg, by = character())

jnds$median <- abs(jnds$l - jnds$bg)
jnds <- jnds[, c(1,3,5)]
colnames(jnds) <- c("abbrevs.x","abbrevs.y","median")
medsLum <- jnds

medsLum$group <- "0"
medsLum[medsLum$abbrevs.x == medsLum$abbrevs.y,"group"] <- "local"
medsLum[medsLum$abbrevs.x != medsLum$abbrevs.y,"group"] <- "foreign"

medsLum <- medsLum %>%
  group_by(abbrevs.y) %>%
  mutate(good_ranks = order(order(median, decreasing=FALSE)))
loc <- medsLum[medsLum$abbrevs.x == medsLum$abbrevs.y,]
png(file = 'output//loc-nonloc//Lum_bc_loc-ranks_centroids.png', width = 1024, height = 768, units = "px")
hist(loc$good_ranks, xlab = "rank", main = "local lithops distance rankings",
     breaks = seq(min(loc$good_ranks), max(loc$good_ranks), length.out = 13))
dev.off()
