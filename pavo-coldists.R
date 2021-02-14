library(pavo)

source("loading-cleaning.R")

unique(df$abbrevs)
ls <- split(df, df$abbrevs)

for(i in 1:length(ls)){
  
  df <- ls[[i]]
  df <- df[,c("abbrevs","roi","mspec","swMean","mwMean","lwMean","lumMean")]
  df <- unite(df,"label",c("abbrevs","mspec","roi"), sep = "__")
  rnames <- df[,1:3]
  rownames(df) <- rnames[,1]
  df <- df[,-1]
  
  # what does 'neural' and 'quantum' noise really refer to?
  dat <- coldist(dat,achromatic = TRUE, subset = c("__a"), qcatch = "Qi", n = c(1,1,1),
                 weber = c(.1,.1,.1), weber.achro = .1, noise = "neural") 
  
  dat <- separate(dat,patch1, into = c("abbrevs.x","mspec.x","roi.x"), sep = "__")
  dat <- separate(dat,patch2, into = c("abbrevs.y","mspec.y","roi.y"), sep = "__")
  dat$substrate.x <- substr(dat$roi.x,0,1)
  dat$substrate.y <- substr(dat$roi.y,0,1)
  jnds <- dat[dat$substrate.x != dat$substrate.y,]
  
  jnds <- jnds %>%  mutate(quartile = ntile(dL, 4))
  closeJND <-jnds[jnds$quartile == 1,]
  
  ls[[i]] <- closeJND %>% group_by(abbrevs.x) %>% summarise(median = median(dS))
}
x <- bind_rows(ls, .id = "abbrevs.x") # why NAs? 
x <- na.omit(x)

png(file = 'output//lum-near-dist-medians_pavo-crossref.png', width = 1024, height = 768, units = "px")
ggplot(x, aes(x = abbrevs.x,y=median)) +
  geom_jitter(data = x, width=0.1,alpha=0.8,color = "red", size =2) + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2,
                                   size = 8, face = "bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")  
dev.off()
