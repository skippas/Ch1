
# BRIGHTNESS 

# counts of data within threshold - very similar to a nearest set of distances

nLength <- df %>% group_by(abbrevs, group) %>% summarise(n = n())
nLength <- nLength[nLength$group == "bg",]

closeJND <- jnds[jnds$dL < 3,]
freqs <- closeJND %>% group_by(abbrevs.y, abbrevs.x,mspec.x, roi.x) %>% summarise(n = n())
freqs <- merge(x = nLength, y = freqs, by.x = "abbrevs", by.y = "abbrevs.y", all.y = TRUE)
freqs$prop.near.plant <- freqs$n.y / freqs$n.x
freqs <- freqs[, !names(freqs) %in% c("group","n.x", "n.y")]
freqs <- freqs %>% group_by(abbrevs, abbrevs.x) %>% summarise(prop.near.pop = mean(prop.near.plant))

freqs$group <- "group"
freqs[freqs$abbrevs.x == freqs$abbrevs,"group"] <- "local"
freqs[freqs$abbrevs.x != freqs$abbrevs,"group"] <- "foreign"
unique(freqs$abbrevs.x)

freqs$abbrevs <- with(freqs, reorder(abbrevs, prop.near.pop, function(x) max(x)))

png(file = 'output//prop-bg_brightness_order-max.png', width = 1024, height = 768, units = "px")
ggplot(freqs, aes(x = abbrevs,y=prop.near.pop)) +
  geom_jitter(data = freqs %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = freqs %>% filter(group == "local"),width=0.1,alpha=0.5,color = "red") + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2, size = 8), axis.title.x = element_blank()) + theme_bw()
dev.off()

# order by local match distance from best overall match
freqs <- freqs %>% 
  group_by(abbrevs) %>% 
  mutate(all.max = max(prop.near.pop))
x <- freqs %>% group_by(abbrevs, abbrevs.x) %>%
  filter(abbrevs == abbrevs.x) %>%
  mutate(local.max = max(prop.near.pop))
x <- x[,c("abbrevs", "local.max")]
freqs <- merge(freqs, x, by = "abbrevs", all.x = TRUE)
freqs$local.diff <- freqs$all.max - freqs$local.max

freqs$abbrevs <- with(freqs, reorder(abbrevs, local.diff, function(x) min(x)))

png(file = 'output//prop-bg_brightness.png', width = 1024, height = 768, units = "px")
ggplot(freqs, aes(x = abbrevs,y=prop.near.pop)) +
  geom_jitter(data = freqs %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = freqs %>% filter(group == "local"),width=0.1,alpha=0.5,color = "red")  + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2, size = 8), axis.title.x = element_blank())
dev.off()


# COLOUR 

# counts of data within threshold - very similar to other things i've done
# But its important to account for range / variance of Lithops by taking the average count, per plant 

# looks like theres duplicate lithops names causing an issue. need to unite with mspec and then make unique - data cleaning section
nLength <- df %>% group_by(abbrevs, group) %>% summarise(n = n())
nLength <- nLength[nLength$group == "bg",]

closeJND <- jnds[jnds$dS < 3,]
freqs <- closeJND %>% group_by(abbrevs.y, abbrevs.x,mspec.x, roi.x) %>% summarise(n = n())
freqs <- merge(x = nLength, y = freqs, by.x = "abbrevs", by.y = "abbrevs.y", all.y = TRUE)
freqs$prop.near.plant <- freqs$n.y / freqs$n.x
freqs <- freqs[, !names(freqs) %in% c("group","n.x", "n.y")]
freqs <- freqs %>% group_by(abbrevs, abbrevs.x) %>% summarise(prop.near.pop = mean(prop.near.plant))

freqs$group <- "group"
freqs[freqs$abbrevs.x == freqs$abbrevs,"group"] <- "local"
freqs[freqs$abbrevs.x != freqs$abbrevs,"group"] <- "foreign"
unique(freqs$abbrevs.x)

freqs$abbrevs <- with(freqs, reorder(abbrevs, prop.near.pop, function(x) max(x)))

png(file = 'output//prop-bg_colour_order-max.png', width = 1024, height = 768, units = "px")
ggplot(freqs, aes(x = abbrevs,y=prop.near.pop)) +
  geom_jitter(data = freqs %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = freqs %>% filter(group == "local"),width=0.1,alpha=0.5,color = "red") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2, size = 8), axis.title.x = element_blank())
dev.off()

# order by local match distance from best overall match
freqs <- freqs %>% 
  group_by(abbrevs) %>% 
  mutate(all.max = max(prop.near.pop))
x <- freqs %>% group_by(abbrevs, abbrevs.x) %>%
  filter(abbrevs == abbrevs.x) %>%
  mutate(local.max = max(prop.near.pop))
x <- x[,c("abbrevs", "local.max")]
freqs <- merge(freqs, x, by = "abbrevs", all.x = TRUE)
freqs$local.diff <- freqs$all.max - freqs$local.max

freqs$abbrevs <- with(freqs, reorder(abbrevs, local.diff, function(x) min(x)))

png(file = 'output//prop-bg_colour.png', width = 1024, height = 768, units = "px")
ggplot(freqs, aes(x = abbrevs,y=prop.near.pop)) +
  geom_jitter(data = freqs %>% filter(group == "foreign"),width=0.1,alpha=0.2,color = "black") + 
  geom_jitter(data = freqs %>% filter(group == "local"),width=0.1,alpha=0.5,color = "red")  + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2, size = 8), axis.title.x = element_blank())
dev.off()






