# https://web.stanford.edu/class/bios221/labs/simulation/Lab_3_simulation.html
library(tidyverse)

scen1 <- rnorm(2000000, 0.5,.03)
scen1 <- data.frame(reflectance = scen1,
                    background_patch = rep("patch 1", 2e6),
                    scenario = "scenario 1")

patch1=rnorm(2000000, 0.75,.03)
patch2=rnorm(1000000, 0.25, .03)
scen2 <- data.frame(reflectance=c(patch1, patch2),
           background_patch=rep(c("patch 1", "patch 2"),
                      c(2e6, 1e6),),
           scenario = "scenario 2")

patch1=rnorm(1000000, 0.5,.03)
patch2=rnorm(1000000, 0.45, .03)
overall = c(patch1, patch2)
scen3 <- data.frame(reflectance=c(overall, patch1, patch2),
                    background_patch=rep(c("overall", "patch 1", "patch 2"),
                                         c(2e6, 1e6, 1e6)),
                    scenario = "scenario 3")

scenarios <- rbind(scen1, scen2, scen3) 
scen_means <- scenarios %>%
  filter(background_patch != "overall") %>%
  group_by(scenario) %>% summarise(mean = mean(reflectance))


png("output\\intro-bgMatching-scenarios.png", type = 'cairo', units = 'mm',
    width = 200, height = 100, res = 300)
ggplot(scenarios,aes(x=reflectance))+ 
  stat_bin(aes(fill = background_patch),
           position="identity",
           binwidth=0.0025, alpha=0.5)+
  xlim(c(0,1))+
  geom_vline(data = scen_means,aes(xintercept = mean), linetype = 2)+
  facet_grid(~scenario)+
  theme_bw()
dev.off()

