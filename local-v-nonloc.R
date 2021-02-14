source("functions//coldist_effic.R")
source("functions//nearest_qtile.R")
# brightness - mean of the nearest quartile of distances local and non local

jnds <- colDistEffic()
medsLum <- nearestQtile(df = jnds,visfeature = "dL", ntiles = 4, nearest_ntile = 1) 

medsLum <- na.omit(medsLum)
medsLum$group <- "0"
medsLum[medsLum$abbrevs.x == medsLum$abbrevs.y,"group"] <- "local"
medsLum[medsLum$abbrevs.x != medsLum$abbrevs.y,"group"] <- "foreign"

# ORDER
medsLum$abbrevs.y <- as.factor(medsLum$abbrevs.y)
ascend <- medsLum[medsLum$group == "local",]
ascend$abbrevs.y <- fct_reorder(ascend$abbrevs.y, ascend$median)
medsLum$abbrevs.y <- fct_relevel(medsLum$abbrevs.y, levels(ascend$abbrevs.y))

png(file = 'output//lum-near-dist-medians.png', width = 1024, height = 768, units = "px")
ggplot(medsLum, aes(x = abbrevs.y,y=median)) +
  geom_jitter(data = medsLum %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = medsLum %>% filter(group == "local"),width=0.1,alpha=0.8,color = "red", size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=0.95,vjust=0.9, size = 8, face = "bold")) + 
  #theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")
dev.off()

## rank distribution
medsLum <- medsLum %>% mutate(ranking = rank(median))
loc <- medsLum[medsLum$abbrevs.x == medsLum$abbrevs.y,]
nonloc <- medsLum[medsLum$abbrevs.x != medsLum$abbrevs.y,]
png(file = 'output//lum-local-dist-ranks.png', width = 1024, height = 768, units = "px")
hist(loc$ranking, xlab = "rank", main = "local lithops distance rankings",
     breaks = seq(min(loc$ranking), max(loc$ranking), length.out = 13))
dev.off()
## distribution of median distances
hist(loc$median, seq(0,70,2))
hist(nonloc$median, seq(0,70,2))
ggplot(medsLum) + 
  geom_freqpoly(aes(x = median, colour = group), stat = "density") +
  xlab("median JND") + geom_vline(xintercept = 1, linetype = "dotted") +
  annotate("text", x=4.0, y=1.1, label=paste("1 JND"), size=4, fontface ="bold") +
  theme_bw()

p1 <- ggplot(data = loc) + geom_histogram(aes(x = median),alpha = 0.8, fill = "black", breaks = seq(0,70,2)) + 
  geom_vline(xintercept = median(loc$median), linetype = "dashed") + 
  theme_bw()
p2 <- ggplot(data = nonloc) + geom_histogram(aes(x = median),alpha = 0.8, fill = "black", breaks = seq(0,70,2)) + 
  geom_vline(xintercept = median(nonloc$median), linetype = "dashed") + 
  theme_bw()
grid.arrange(p1,p2)

# COLOUR 

medsCol <- nearestQtile(visfeature = "dS", ntiles = 4, nearest_ntile = 1) 

medsCol <- na.omit(medsCol)
medsCol$group <- "0"
medsCol[medsCol$abbrevs.x == medsCol$abbrevs.y,"group"] <- "local"
medsCol[medsCol$abbrevs.x != medsCol$abbrevs.y,"group"] <- "foreign"

# ORDER
medsCol$abbrevs.y <- as.factor(medsCol$abbrevs.y)
ascend <- medsCol[medsCol$group == "local",]
ascend$abbrevs.y <- fct_reorder(ascend$abbrevs.y, ascend$median)
medsCol$abbrevs.y <- fct_relevel(medsCol$abbrevs.y, levels(ascend$abbrevs.y))

png(file = 'output//color-near-dist-medians.png', width = 1024, height = 768, units = "px")
ggplot(medsCol, aes(x = abbrevs.y,y=median)) +
  geom_jitter(data = medsCol %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + theme_bw() +
  geom_jitter(data = medsCol %>% filter(group == "local"),width=0.1,alpha=0.8,color = "red", size =2) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 12, face = "bold"),
        axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")  
dev.off()


## rank distribution
medsCol <- medsCol %>% mutate(ranking = rank(median))
loc <- medsCol[medsCol$abbrevs.x == medsCol$abbrevs.y,]
nonloc <- medsCol[medsCol$abbrevs.x != medsCol$abbrevs.y,]
png(file = 'output//col-local-dist-ranks.png', width = 1024, height = 768, units = "px")
hist(loc$ranking, xlab = "rank", main = "local lithops distance rankings",
     breaks = seq(min(loc$ranking), max(loc$ranking), length.out = 13))
dev.off()
## distribution of median distances
hist(loc$median, seq(0,10,1))
hist(nonloc$median, seq(0,20,1))
ggplot(medsCol) + 
  geom_freqpoly(aes(x = median, colour = group), stat = "density") +
  xlab("median JND") + geom_vline(xintercept = 1, linetype = "dotted") +
  annotate("text", x=1.7, y=1.1, label=paste("1 JND"), size=4, fontface ="bold") +
  theme_bw()

p1 <- ggplot(data = loc) + geom_histogram(aes(x = median),alpha = 0.8, fill = "black", breaks = seq(0,70,2)) + 
  geom_vline(xintercept = median(loc$median), linetype = "dashed") + 
  theme_bw()
p2 <- ggplot(data = nonloc) + geom_histogram(aes(x = median),alpha = 0.8, fill = "black", breaks = seq(0,70,2)) + 
  geom_vline(xintercept = median(nonloc$median), linetype = "dashed") + 
  theme_bw()
grid.arrange(p1,p2)

# unlikely local combinations - worst foreign combinations
poorMatch <- medsCol[medsCol$median > 3,]
poorMatch <- medsLum[medsLum$median > 8, ]

