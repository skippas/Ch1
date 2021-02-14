# brightness
rm(list = ls())
source("functions//coldist_effic.R")
source("functions//nearest_qtile.R")

jnds <- colDistEffic(c("b","c"))
medsLumBC <- nearestQtile(df = jnds,visfeature = "dL", ntiles = 4, nearest_ntile = 4) 
jnds <- colDistEffic("b")
medsLumB <- nearestQtile(df = jnds,visfeature = "dL", ntiles = 4, nearest_ntile = 4) 
jnds <- colDistEffic("c")
medsLumC <- nearestQtile(df = jnds,visfeature = "dL", ntiles = 4, nearest_ntile = 4) 
colnames(medsLumB)[which(names(medsLumB) == "median")] <- "median_b"
colnames(medsLumC)[which(names(medsLumC) == "median")] <- "median_c"
colnames(medsLumBC)[which(names(medsLumBC) == "median")] <- "median_bc"
medsLum <- merge(medsLumB[,c("abbrevs.x","abbrevs.y","median_b")], medsLumC[,c("abbrevs.x","abbrevs.y","median_c")])
medsLum <- merge(medsLum[,c("abbrevs.x","abbrevs.y","median_b","median_c")], medsLumBC[,c("abbrevs.x","abbrevs.y","median_bc")])
medsLum$group <- "foreign"
medsLum[medsLum$abbrevs.x == medsLum$abbrevs.y,"group"] <- "local"
medsLum$abbrevs.y <- as.factor(medsLum$abbrevs.y)
medsLum <- pivot_longer(medsLum, cols = c("median_b", "median_c", "median_bc"),
                        names_to = "substrate", values_to = "median")

ascend <- medsLum %>% filter(group == "local", substrate == "median_b") %>%
  group_by(abbrevs.x, abbrevs.y) %>% slice(which.min(median))
ascend$abbrevs.y <- fct_reorder(ascend$abbrevs.y, ascend$median)
medsLum$abbrevs.y <- fct_relevel(medsLum$abbrevs.y, levels(ascend$abbrevs.y))

png(file = "output//loc-nonloc//LUM_b.png", width = 1024, height = 768, units = "px")
ggplot(medsLum, aes(x = abbrevs.y,y=median)) +
  geom_jitter(data = medsLum %>% filter(group == "foreign", substrate == "median_b"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = medsLum %>% filter(group == "local", substrate == "median_b"),width=0.1,alpha=0.8,color = "red", size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=0.95,vjust=0.9, size = 8, face = "bold")) + 
  #theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")
dev.off()

png(file = "output//loc-nonloc//LUM_b_c.png", width = 1024, height = 768, units = "px")
ggplot(medsLum, aes(x = abbrevs.y,y=median)) +
  geom_jitter(data = medsLum %>% filter(group == "foreign", substrate == "median_b"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = medsLum %>% filter(group == "local", substrate == "median_b"),width=0.1,alpha=0.8,color = "red", size = 2) +
  geom_jitter(data = medsLum %>% filter(group == "local", substrate == "median_c"),width=0.1,alpha=0.8,color = "blue", size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=0.95,vjust=0.9, size = 8, face = "bold")) + 
  #theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")
dev.off()

minVal <- medsLum %>% filter(substrate == "median_bc") %>%
  group_by(abbrevs.y) %>% slice(which.min(median))
locVal <- medsLum %>% filter(group == "local", substrate == "median_bc") 
ascend <- rbind(locVal, minVal)
ascend <- aggregate(median~abbrevs.y, ascend, FUN= function(i)max(i) - min(i)) 
ascend$abbrevs.y <- fct_reorder(ascend$abbrevs.y, ascend$median)
medsLum$abbrevs.y <- fct_relevel(medsLum$abbrevs.y, levels(ascend$abbrevs.y))

png(file = "output//loc-nonloc//LUM_bc.png", width = 1024, height = 768, units = "px")
ggplot(medsLum, aes(x = abbrevs.y,y=median)) +
  geom_jitter(data = medsLum %>% filter(group == "foreign", substrate == "median_bc"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = medsLum %>% filter(group == "local", substrate == "median_bc"),width=0.1,alpha=0.8,color = "red", size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=0.95,vjust=0.9, size = 8, face = "bold")) + 
  #theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")
dev.off()

# RANKS
ranks <- medsLum %>% filter(substrate == "median_bc") %>% 
  group_by(abbrevs.y) %>% mutate(ranking = rank(median))
locRanks <- ranks[ranks$abbrevs.x == ranks$abbrevs.y,]
png(file = 'output//loc-nonloc//Lum_bc_loc-ranks.png', width = 1024, height = 768, units = "px")
hist(locRanks$ranking, xlab = "rank", main = "local lithops distance rankings",
     breaks = seq(min(locRanks$ranking), max(locRanks$ranking), length.out = 13))
dev.off()

ranks <- medsLum %>% filter(substrate == "median_b") %>% 
  group_by(abbrevs.y) %>% mutate(ranking = rank(median))
locRanks <- ranks[ranks$abbrevs.x == ranks$abbrevs.y,]
png(file = 'output//loc-nonloc//Lum_b_loc-ranks.png', width = 1024, height = 768, units = "px")
hist(locRanks$ranking, xlab = "rank", main = "local lithops distance rankings",
     breaks = seq(min(locRanks$ranking), max(locRanks$ranking), length.out = 13))
dev.off()


# COLOUR
rm(list = ls())
source("functions//coldist_effic.R")
source("functions//nearest_qtile.R")

jnds <- colDistEffic(c("b","c"))
medsColBC <- nearestQtile(df = jnds,visfeature = "dS", ntiles = 4, nearest_ntile = 4) 
jnds <- colDistEffic("b")
medsColB <- nearestQtile(df = jnds,visfeature = "dS", ntiles = 4, nearest_ntile = 4) 
jnds <- colDistEffic("c")
medsColC <- nearestQtile(df = jnds,visfeature = "dS", ntiles = 4, nearest_ntile = 4) 
colnames(medsColB)[which(names(medsColB) == "median")] <- "median_b"
colnames(medsColC)[which(names(medsColC) == "median")] <- "median_c"
colnames(medsColBC)[which(names(medsColBC) == "median")] <- "median_bc"
medsCol <- merge(medsColB[,c("abbrevs.x","abbrevs.y","median_b")], medsColC[,c("abbrevs.x","abbrevs.y","median_c")])
medsCol <- merge(medsCol[,c("abbrevs.x","abbrevs.y","median_b","median_c")], medsColBC[,c("abbrevs.x","abbrevs.y","median_bc")])
medsCol$group <- "foreign"
medsCol[medsCol$abbrevs.x == medsCol$abbrevs.y,"group"] <- "local"
medsCol$abbrevs.y <- as.factor(medsCol$abbrevs.y)
medsCol <- pivot_longer(medsCol, cols = c("median_b", "median_c", "median_bc"),
                        names_to = "substrate", values_to = "median")

ascend <- medsCol %>% filter(group == "local", substrate == "median_b") %>%
  group_by(abbrevs.x, abbrevs.y) %>% slice(which.min(median))
ascend$abbrevs.y <- fct_reorder(ascend$abbrevs.y, ascend$median)
medsCol$abbrevs.y <- fct_relevel(medsCol$abbrevs.y, levels(ascend$abbrevs.y))

png(file = "output//loc-nonloc//Col_b_c.png", width = 1024, height = 768, units = "px")
ggplot(medsCol, aes(x = abbrevs.y,y=median)) +
  geom_jitter(data = medsCol %>% filter(group == "foreign", substrate == "median_b"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = medsCol %>% filter(group == "local", substrate == "median_b"),width=0.1,alpha=0.8,color = "red", size = 2) +
  geom_jitter(data = medsCol %>% filter(group == "local", substrate == "median_c"),width=0.1,alpha=0.8,color = "blue", size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=0.95,vjust=0.9, size = 8, face = "bold")) + 
  #theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")
dev.off()

minVal <- medsCol %>% filter(substrate == "median_bc") %>%
  group_by(abbrevs.y) %>% slice(which.min(median))
locVal <- medsCol %>% filter(group == "local", substrate == "median_bc") 
ascend <- rbind(locVal, minVal)
ascend <- aggregate(median~abbrevs.y, ascend, FUN= function(i)max(i) - min(i)) 
ascend$abbrevs.y <- fct_reorder(ascend$abbrevs.y, ascend$median)
medsCol$abbrevs.y <- fct_relevel(medsCol$abbrevs.y, levels(ascend$abbrevs.y))

png(file = "output//loc-nonloc//Col_bc.png", width = 1024, height = 768, units = "px")
ggplot(medsCol, aes(x = abbrevs.y,y=median)) +
  geom_jitter(data = medsCol %>% filter(group == "foreign", substrate == "median_bc"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = medsCol %>% filter(group == "local", substrate == "median_bc"),width=0.1,alpha=0.8,color = "red", size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=0.95,vjust=0.9, size = 8, face = "bold")) + 
  #theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")
dev.off()

# RANKS
ranks <- medsCol %>% filter(substrate == "median_bc") %>% 
  group_by(abbrevs.y) %>% mutate(ranking = rank(median))
locRanks <- ranks[ranks$abbrevs.x == ranks$abbrevs.y,]
png(file = 'output//loc-nonloc//Col_bc_loc-ranks.png', width = 1024, height = 768, units = "px")
hist(locRanks$ranking, xlab = "rank", main = "local lithops distance rankings",
     breaks = seq(min(locRanks$ranking), max(locRanks$ranking), length.out = 13))
dev.off()


locVnonloc <- function(sub = c("c","b"), visfeat = "dL") {
  jnds <- colDistEffic(sub)
  medsLum <- nearestQtile(data = jnds,visfeature = visfeat, ntiles = 4, nearest_ntile = 4) 
  medsLum$group <- "foreign"
  medsLum[medsLum$abbrevs.x == medsLum$abbrevs.y,"group"] <- "local"

  filename <- c("output//",visfeat,"-near-dist-medians.png")
  filename <- paste0(filename, collapse = "")
  png(file = filename, width = 1024, height = 768, units = "px")
  ggplot(medsLum, aes(x = abbrevs.y,y=median)) +
    geom_jitter(data = medsLum %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + 
    geom_jitter(data = medsLum %>% filter(group == "local"),width=0.1,alpha=0.8,color = "red", size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=0.95,vjust=0.9, size = 8, face = "bold")) + 
    #theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
    labs(y = "lithops median distance to background", x = "background populations")
  dev.off()
  
}
x <- locVnonloc(sub = "c")

jnds <- colDistEffic("b")
medsLum_b <- nearestQtile(df = jnds,visfeature = "dL", ntiles = 4, nearest_ntile = 4) 
medsLum_b$group <- "foreign"
medsLum_b[medsLum_b$abbrevs.x == medsLum_b$abbrevs.y,"group"] <- "local"
colnames(medsLum_b)[which(names(medsLum_b) == "median")] <- "median_b"
locMedsLum_b <- medsLum_b[medsLum_b$abbrevs.x == medsLum_b$abbrevs.y,]
frnMedsLum_b <- medsLum_b[medsLum_b$abbrevs.x != medsLum_b$abbrevs.y,] 

jnds <- colDistEffic("c")
medsLum_c <- nearestQtile(df = jnds,visfeature = "dL", ntiles = 4, nearest_ntile = 4) 
medsLum_c$group <- "foreign"
medsLum_c[medsLum_c$abbrevs.x == medsLum_c$abbrevs.y,"group"] <- "local"
colnames(medsLum_b)[which(names(medsLum_b) == "median")] <- "median_c"
locMedsLum_c <- medsLum_c[medsLum_c$abbrevs.x == medsLum_c$abbrevs.y,] 
frnMedsLum_c <- medsLum_c[medsLum_c$abbrevs.x != medsLum_c$abbrevs.y,] 

frnMedsLum <- merge(frnMedsLum_c[,c("abbrevs.x", "abbrevs.y", "median_c")],
                    frnMedsLum_b[,c("abbrevs.x", "abbrevs.y", "median_b")],
                    by.x = c("abbrevs.x", "abbrevs.y"), by.y = c("abbrevs.x", "abbrevs.y"))
frnMedsLum <- pivot_longer(frnMedsLum, cols = c("median_c", "median_b"),
                           names_to = "substrate", values_to = "median")
#frnMedsLum <- frnMedsLum %>% group_by(abbrevs.x,abbrevs.y) %>% slice(which.min(median))
frnMedsLum$group <- "foreign"

locMedsLum <- merge(locMedsLum_c[,c("abbrevs.x", "abbrevs.y", "median_c")],
                    locMedsLum_b[,c("abbrevs.x", "abbrevs.y", "median_b")],
                    by.x = c("abbrevs.x", "abbrevs.y"), by.y = c("abbrevs.x", "abbrevs.y"))
locMedsLum <- pivot_longer(locMedsLum, cols = c("median_c", "median_b"),
                           names_to = "substrate", values_to = "median")
locMedsLum$group <- "local"

bestSub <- locMedsLum %>% group_by(abbrevs.x, abbrevs.y) %>% slice(which.min(median))
bestSub$id <- paste(bestSub$abbrevs.y, bestSub$substrate)
frnMedsLum$id <- paste(frnMedsLum$abbrevs.y, frnMedsLum$substrate)
frnMedsLum <- frnMedsLum[frnMedsLum$id %in% bestSub$id, , drop = FALSE]
frnMedsLum <- subset(frnMedsLum ,select = -c(id))

medsLum <- rbind(locMedsLum, frnMedsLum)

# ORDER
locMedsLum$abbrevs.y <- as.factor(locMedsLum$abbrevs.y)
ascend <- locMedsLum %>% group_by(abbrevs.x, abbrevs.y) %>% slice(which.min(median))
ascend$abbrevs.y <- fct_reorder(ascend$abbrevs.y, ascend$median)
medsLum$abbrevs.y <- fct_relevel(medsLum$abbrevs.y, levels(ascend$abbrevs.y))

png(file = 'output//local-nonloc-meds_lum_by-sub.png', width = 1024, height = 768, units = "px")
ggplot() +
  geom_jitter(data = medsLum %>% filter(group == "foreign"), aes(x = abbrevs.y, y = median), width=0.1,alpha=0.2,color = "black") +
  geom_jitter(data = medsLum %>% filter(group == "local" & substrate == "median_b"), aes(x = abbrevs.y, y = median), width=0.1,alpha=0.8,color = "red", size = 2) +
  geom_jitter(data = medsLum %>% filter(group == "local" & substrate == "median_c"), aes(x = abbrevs.y, y = median), width=0.1,alpha=0.8,color = "red", size = 2) +
  theme_bw() + geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 12, face = "bold"),
  #axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")  
dev.off()

## rank distribution
temp <- medsLum %>% filter(group == "local") %>% group_by(abbrevs.x, abbrevs.y) %>%
  slice(which.min(median))
medsLum <- medsLum %>% filter(group == "foreign")
medsLum <- rbind(medsLum, temp)

medsLum <- medsLum %>% group_by(abbrevs.y) %>% mutate(ranking = rank(median))
loc <- medsLum[medsLum$abbrevs.x == medsLum$abbrevs.y,]
nonloc <- medsLum[medsLum$abbrevs.x != medsLum$abbrevs.y,]
png(file = 'output//lum-local-dist-ranks.png', width = 1024, height = 768, units = "px")
hist(loc$ranking, xlab = "rank", main = "local lithops distance rankings",
     breaks = seq(min(loc$ranking), max(loc$ranking), length.out = 13))
dev.off()

# COLOUR 

jnds <- colDistEffic("b")
medsCol_b <- nearestQtile(df = jnds,visfeature = "dS", ntiles = 4, nearest_ntile = 4) 

medsCol_b <- na.omit(medsCol_b)
medsCol_b$group <- "0"
medsCol_b[medsCol_b$abbrevs.x == medsCol_b$abbrevs.y,"group"] <- "local"
medsCol_b[medsCol_b$abbrevs.x != medsCol_b$abbrevs.y,"group"] <- "foreign"
locMedsCol_b <- medsCol_b[medsCol_b$abbrevs.x == medsCol_b$abbrevs.y,]
frnMedsCol_b <- medsCol_b[medsCol_b$abbrevs.x != medsCol_b$abbrevs.y,] 

jnds <- colDistEffic("c")
medsCol_c <- nearestQtile(df = jnds,visfeature = "dS", ntiles = 4, nearest_ntile = 4) 

medsCol_c <- na.omit(medsCol_c)
medsCol_c$group <- "0"
medsCol_c[medsCol_c$abbrevs.x == medsCol_c$abbrevs.y,"group"] <- "local"
medsCol_c[medsCol_c$abbrevs.x != medsCol_c$abbrevs.y,"group"] <- "foreign"
locMedsCol_c <- medsCol_c[medsCol_c$abbrevs.x == medsCol_c$abbrevs.y,] 
frnMedsCol_c <- medsCol_c[medsCol_c$abbrevs.x != medsCol_c$abbrevs.y,] 

frnMedsCol_b$median_b <- frnMedsCol_b$median
frnMedsCol_c$median_c <- frnMedsCol_c$median
frnMedsCol <- merge(frnMedsCol_c[,c("abbrevs.x", "abbrevs.y", "median_c")],
                    frnMedsCol_b[,c("abbrevs.x", "abbrevs.y", "median_b")],
                    by.x = c("abbrevs.x", "abbrevs.y"), by.y = c("abbrevs.x", "abbrevs.y"))
frnMedsCol <- pivot_longer(frnMedsCol, cols = c("median_c", "median_b"),
                           names_to = "substrate", values_to = "median")
#frnMedsLum <- frnMedsLum %>% group_by(abbrevs.x,abbrevs.y) %>% slice(which.min(median))
frnMedsCol$group <- "foreign"

locMedsCol_b$median_b <- locMedsCol_b$median
locMedsCol_c$median_c <- locMedsCol_c$median
locMedsCol <- merge(locMedsCol_c[,c("abbrevs.x", "abbrevs.y", "median_c")],
                    locMedsCol_b[,c("abbrevs.x", "abbrevs.y", "median_b")],
                    by.x = c("abbrevs.x", "abbrevs.y"), by.y = c("abbrevs.x", "abbrevs.y"))
locMedsCol <- pivot_longer(locMedsCol, cols = c("median_c", "median_b"),
                           names_to = "substrate", values_to = "median")
locMedsCol$group <- "local"

bestSub <- locMedsCol %>% group_by(abbrevs.x, abbrevs.y) %>% slice(which.min(median))
bestSub$id <- paste(bestSub$abbrevs.y, bestSub$substrate)
frnMedsCol$id <- paste(frnMedsCol$abbrevs.y, frnMedsCol$substrate)
frnMedsCol <- frnMedsCol[frnMedsCol$id %in% bestSub$id, , drop = FALSE]
frnMedsCol <- subset(frnMedsCol ,select = -c(id))

medsCol <- rbind(locMedsCol, frnMedsCol)

# ORDER
locMedsCol$abbrevs.y <- as.factor(locMedsCol$abbrevs.y)
ascend <- locMedsCol %>% group_by(abbrevs.x, abbrevs.y) %>% slice(which.min(median))
ascend$abbrevs.y <- fct_reorder(ascend$abbrevs.y, ascend$median)
medsCol$abbrevs.y <- fct_relevel(medsCol$abbrevs.y, levels(ascend$abbrevs.y))

png(file = 'output//local-nonloc//col_by-sub.png', width = 1024, height = 768, units = "px")
ggplot() +
  geom_jitter(data = medsCol %>% filter(group == "foreign"), aes(x = abbrevs.y, y = median), width=0.1,alpha=0.2,color = "black") +
  geom_jitter(data = medsCol %>% filter(group == "local" & substrate == "median_b"), aes(x = abbrevs.y, y = median), width=0.1,alpha=0.8,color = "red", size = 2) +
  geom_jitter(data = medsCol %>% filter(group == "local" & substrate == "median_c"), aes(x = abbrevs.y, y = median), width=0.1,alpha=0.8,color = "blue", size = 2) +
  theme_bw() + geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 12, face = "bold"),
        #axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")  
dev.off()

## rank distribution
temp <- medsCol %>% filter(group == "local") %>% group_by(abbrevs.x, abbrevs.y) %>%
  slice(which.min(median))
medsCol <- medsCol %>% filter(group == "foreign")
medsCol <- rbind(medsCol, temp)

medsCol <- medsCol %>% group_by(abbrevs.y) %>% mutate(ranking = rank(median))
loc <- medsCol[medsCol$abbrevs.x == medsCol$abbrevs.y,]
nonloc <- medsCol[medsCol$abbrevs.x != medsCol$abbrevs.y,]
png(file = 'output//col-local-dist-ranks.png', width = 1024, height = 768, units = "px")
hist(loc$ranking, xlab = "rank", main = "local lithops distance rankings",
     breaks = seq(min(loc$ranking), max(loc$ranking), length.out = 13))
dev.off()

## rank distribution
x <- near.meds %>% mutate(ranking = rank(median))
hist(x$ranking[x$abbrevs == x$abbrevs.y])
x <- x[x$abbrevs.x == x$abbrevs.y,]
png(file = 'output//lum-local-dist-ranks.png', width = 1024, height = 768, units = "px")
hist(x$ranking, xlab = "rank", main = "local lithops distance rankings",
     breaks = seq(min(x$ranking), max(x$ranking), length.out = 13))
dev.off()

