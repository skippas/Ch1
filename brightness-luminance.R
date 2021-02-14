library(ggrepel)
library(gridExtra)
source("Loading-cleaning.R")

df <- df[df$lumMean > 0,] 
# scale luminance. Not 100 percent about this, looks like
# what is done in pavo
df$lumMean <- df$lumMean*10000 
dat <- df %>% group_by(abbrevs, substrate) %>%
  summarise(meanLum = mean(lumMean), logLum = mean(log(lumMean)))
dat <- pivot_wider(dat, names_from = substrate, values_from = c(logLum, meanLum))

p1 <- ggplot(data = dat, aes(x = meanLum_b, y = meanLum_a)) +
  geom_point(aes(x = meanLum_b, y = meanLum_a)) +
  scale_color_gradient() +
  geom_abline(slope = 1, lty = 2) + 
  geom_smooth(method='lm', formula= y~x) +
  theme_bw()
p2 <- ggplot(data = dat,aes(x = logLum_b, y = logLum_a)) +
  geom_point(aes(x = logLum_b, y = logLum_a)) +
  scale_color_gradient() +
  geom_abline(slope = 1, lty = 2) + 
  geom_smooth(method='lm', formula= y~x) +
  theme_bw()
p3 <- ggplot(data = dat,aes(x = logLum_c, y = logLum_a)) +
  geom_point(aes(x = logLum_c, y = logLum_a)) +
  scale_color_gradient() +
  geom_abline(slope = 1, lty = 2) + 
  geom_smooth(method='lm', formula= y~x) +
  theme_bw()

png(file = 'output//Qcatch-logQcatch.png', width = 1024, height = 768, units = "px")
grid.arrange(p1, p2, ncol=2)
dev.off()

mod1 <- lm(dat$logLum_a~dat$logLum_b)
mod2 <- lm(dat$logLum_a~dat$logLum_c)
mod3 <- lm(dat$logLum_a~dat$logLum_b + dat$logLum_c+ dat$logLum_b*dat$logLum_c)

# by group 

dat.bg <- df %>% group_by(abbrevs, group) %>%
  summarise(meanLum = mean(lumMean), logLum = mean(log(lumMean)))
dat.bg <- pivot_wider(dat.bg, names_from = group, values_from = c(logLum, meanLum))

p4 <- ggplot(data = dat.bg, aes(x = logLum_bg, y = logLum_lithops)) +
  geom_point(aes(x = logLum_bg, y = logLum_lithops)) +
  scale_color_gradient() +
  geom_abline(slope = 1, lty = 2) + 
  geom_smooth(method='lm', formula= y~x) +
  theme_bw()
png(file = 'output//lum_reg_l-r-s.png', width = 1024, height = 768, units = "px")
grid.arrange(p2, p3,p4, ncol=3)
dev.off()

mod4 <- lm(dat$logLum_lithops~dat$logLum_bg)
summary(mod1)
summary(mod2)
summary(mod3)

# residual explanations...
# matching soil instead
dat <- df %>% group_by(abbrevs, substrate) %>%
  summarise(meanLum = mean(lumMean), logLum = mean(log(lumMean)))

resids <- names(sort(mod1$residuals)[1:15])
resids <- dat[dat$substrate ==  "a", "abbrevs"][resids,]
resids <- dat[dat$abbrevs %in% resids$abbrevs,]
resids <- resids[, -which(names(resids) %in% c("group", "meanLum"))]

sub.prop <- df %>% count(abbrevs, substrate)
sub.prop <- sub.prop %>% pivot_wider(names_from = substrate, values_from = n)
sub.prop$rock_cover <- sub.prop$b / (sub.prop$b + sub.prop$c)
sub.prop$soil_cover <- sub.prop$c / (sub.prop$b + sub.prop$c)
sub.prop <- arrange(sub.prop, as.factor(abbrevs)) 
sub.prop <- sub.prop[, c("abbrevs","soil_cover")]

resids <- merge(resids, sub.prop, by = "abbrevs", all.x = TRUE)
resids <- pivot_wider(resids, names_from = substrate, values_from = logLum)

png(file = 'output//lum_rock-resids-soil-cover.png', width = 1024, height = 768, units = "px")
ggplot(data = resids) +
  geom_point(aes(x = c, y = a, col = soil_cover)) +
  scale_color_gradient() +
  geom_text_repel(aes(x = c, y = a, label = abbrevs)) +
  geom_abline(slope = 1) + theme_bw()
dev.off()

# variability ?
b.dom <- sub.prop[sub.prop$soil_cover < 0.5, "abbrevs"]
vary <- df %>% filter(abbrevs %in% b.dom$abbrevs, substrate == "b") %>%
  group_by(abbrevs) %>% summarise(b.vary = var(log(lumMean)))

dat <- dat %>% ungroup() %>% select(-meanLum) %>% 
  pivot_wider(names_from = substrate, values_from = logLum)
vary <- merge(dat,vary, by = "abbrevs")

png(file = 'output//lum_rock-resids-variab.png', width = 1024, height = 768, units = "px")
ggplot(data = vary) +
  geom_point(aes(x = b, y = a, col = b.vary)) +
  scale_color_gradient() +
  #geom_text_repel(aes(x = c, y = a, label = abbrevs)) +
  geom_abline(slope = 1) + theme_bw()
dev.off()

# soil explanation?
vary <- merge(vary, sub.prop[, c("abbrevs", "soil_cover"), by = "abbrevs"])
ggplot(data = vary) +
  geom_point(aes(x = c, y = a, col = soil_cover)) +
  scale_color_gradient() +
  geom_text_repel(aes(x = c, y = a, label = abbrevs),fontface = "bold", size = 2) +
  geom_abline(slope = 1) + theme_bw()

# largest residuals w soil cover < 25 percent ? whats going on? Variab?

# Nearest distances / amounts matched

# rock-soil-autocorrelation
dat <- merge(dat, sub.prop[, c("abbrevs", "soil_cover"), by = "abbrevs"])
png("output//rock-soil-autocorrelation", width = 1024, height = 768, units = "px")
ggplot(data = dat) +
  geom_point(aes(x = b, y = c, col = soil_cover)) +
  scale_color_gradient() +
  #geom_text_repel(aes(x = b, y = c, label = abbrevs),fontface = "bold", size = 2) +
  geom_abline(slope = 1) + theme_bw()
dev.off()


# rock-regression-with-labels
dat <- merge(dat, sub.prop[, c("abbrevs", "soil_cover"), by = "abbrevs"])
png("output//rock-regression", width = 1024, height = 768, units = "px")
ggplot(data = dat) +
  geom_point(aes(x = b, y = a, col = soil_cover)) +
  scale_color_gradient() +
  geom_text_repel(aes(x = b, y = a, label = abbrevs),fontface = "bold", size = 2) +
  geom_abline(slope = 1) + theme_bw()
dev.off()



