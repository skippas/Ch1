lith_bright <- df[df$substrate == "a",]
ggplot(lith_bright, aes(pop_spp,lumMean)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


