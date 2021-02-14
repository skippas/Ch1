# run necessary ordering first

# brightness - mean of the nearest quartile of distances local and non local
png(file = 'output//lum-near-dist-medians.png', width = 1024, height = 768, units = "px")
ggplot(near.meds, aes(x = abbrevs.x,y=median)) +
  geom_jitter(data = near.meds %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = near.meds %>% filter(group == "local"),width=0.1,alpha=0.8,color = "red", size = 2) + theme_bw() +
  theme(axis.text.x = element_blank(), axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")
dev.off()

png(file = 'output//lum-near-dist-medians.png', width = 1024, height = 768, units = "px")
ggplot(near.meds, aes(x = abbrevs.x,y=median)) +
  geom_jitter(data = near.meds %>% filter(group == "local"),width=0.1,alpha=0.8,color = "red", size = 2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2, size = 8, face = "bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")
dev.off()

# color

png(file = 'output//color-near-dist-medians.png', width = 1024, height = 768, units = "px")
ggplot(near.meds, aes(x = abbrevs.x,y=median)) +
  geom_jitter(data = near.meds %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + theme_bw() +
  geom_jitter(data = near.meds %>% filter(group == "local"),width=0.1,alpha=0.8,color = "red", size =2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 12, face = "bold"),
        axis.title=element_text(size=20,face="bold")) + 
  labs(y = "lithops median distance to background", x = "background populations")  
dev.off()

# brightness median distance to rock versus soil 
texts <- near.meds.bc %>% filter(median.b > 1 | median.c > 5)
png(file = 'output//brightness-match-rock-v-soil.png', width = 1024, height = 768, units = "px")
ggplot(near.meds.bc, aes(x = median.b, y = median.c)) + geom_point() +
  #geom_text(data = texts, aes(label = abbrevs.x,y = median.c, x = median.b)) +
  theme_bw() + coord_fixed() + 
  labs(y = "Lithops median distances to soil", x = "lithops median distances to rock") +
  theme(axis.title = element_text(size = 20)) 
dev.off()

# colour median distance to rock versus soil 
texts <- near.meds.bc %>% filter(median.b > 1 | median.c > 5)
png(file = 'output//color-match-rock-v-soil.png', width = 1024, height = 768, units = "px")
ggplot(near.meds.bc, aes(x = median.b, y = median.c)) + geom_point() +
  #geom_text(data = texts, aes(label = abbrevs.x,y = median.c, x = median.b)) +
  theme_bw() + coord_fixed(ratio = 1) + 
  labs(y = "Lithops median distances to soil", x = "lithops median distances to rock") +
  theme(axis.title = element_text(size = 20)) 
dev.off()


png(file = 'output//prop-bg_colour.png', width = 1024, height = 768, units = "px")
ggplot(freqs, aes(x = abbrevs,y=prop.near.pop)) +
  geom_jitter(data = freqs %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = freqs %>% filter(group == "local"),width=0.3,alpha=0.5,color = "red") + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2, size = 8), axis.title.x = element_blank()) +
  theme_bw()
dev.off()


