# Filtering large and small sized rocks out

aggregate(df[df$substrate=="a","area"]~df[df$substrate=="a","pop_spp"], FUN=min)

for(i in 1:length(unique(df$pop_spp))){
  lithops_sizes <- df[df$substrate == "a",]
}

lithops_sizes <- lithops_sizes[lithops_sizes$area < 5000,]
ggplot(lithops_sizes, aes(pop_spp,area)) + geom_boxplot()

df <- df[df$area < 5000,]
boxplot(df$area ~ df$substrate)
boxplot()