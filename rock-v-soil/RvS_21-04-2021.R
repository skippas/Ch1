# Calculating and visualizing the distances between rock and soil for
# 1. a nearest subset of luminance distances
# 2. all colour distances
source("nearest-distance-versions.R")
library(ggrepel)
plotF <- function(df,label_df, labelVar, xvar, yvar, ylab, xlab){
  ggplot(df, aes(x = !! sym(xvar), y = !! sym(yvar), label = !! sym(labelVar))) + 
    geom_point() + geom_text_repel(data = label_df) +
    geom_abline(slope = 1, intercept = 0) + 
    theme_bw() + coord_fixed(xlim = c(0,15), ylim = c(0,15)) + 
    labs(y = ylab, x = xlab) +
    theme(axis.title = element_text(size = 10)) 
}
plotList <- vector(mode = "list", length = 2)
labels <- nrGmeanDist[0,]
labels$lithpop <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1',
                               labels$lithpop), ' ')
plotList[[1]] <- plotF(nrGmeanDist, labels, labelVar = "lithpop", xvar = "dL_a_c", 
                  yvar = "dL_a_b", ylab = "rock lum dist", xlab = "soil lum dist")
labels <- nrGmeanDistRelaxed[0,]
plotList[[2]] <- plotF(nrGmeanDistRelaxed, labels, labelVar = "lithpop", xvar = "dL_a_c", 
                     yvar = "dL_a_b", ylab = "rock lum dist", xlab = "soil lum dist")
grid.arrange(plotList[[1]], plotList[[2]], nrow = 1)

png(file = 'output//rock-v-soil//dL-RvS-nrst.png',
    width = 1024, height = 768, units = "px")
plotF(df = nrGmeanDist,labels, labelVar = "lithpop", xvar = "dL_a_c", yvar = "dL_a_b", ylab = "dL rock", xlab = "dL soil")
dev.off()

# 2. all colour distances
allGmeanDists <- gMeanDists(df = df, combine_bg = F,substrates = c("b","c"))
allGmeanDists <- allGmeanDists[allGmeanDists$patch1 == "a" | patch2 == "a"]
allGmeanDists <- pivot_wider(allGmeanDists, values_from = c("dL","dS"),
                          names_from = c("patch1", "patch2"))
allGmeanDists <- allGmeanDists[allGmeanDists$comparison == "local",]

png(file = 'output//rock-v-soil//dS-RvS-all.png',
    width = 1024, height = 768, units = "px")
plotF(df = allGmeanDists, xvar = "dL_a_c", yvar = "dL_a_b", ylab = "dL rock", xlab = "dL soil")
dev.off()

