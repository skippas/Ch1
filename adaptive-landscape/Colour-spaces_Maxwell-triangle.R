library(pavo)
library(tidyverse)
library(gridExtra)
library(grid)
source("loading-cleaning.R")

ls <- split(df, df$pop_spp)
ls <- lapply(ls, function(x) split(x, x$substrate))
ls <- lapply(ls, lapply, function(x) x %>% 
               select(swMean, mwMean, lwMean, lumMean))
ls <- lapply(ls, lapply, colspace, space = "tri")
#lapply(names(ls), function(x) plot(ls[[x]], pch = 21, cex = 1, main = x))

plots <- list(1:59)
pdf("OUTPUT//colour-space_maxwell-triangle_02.pdf")
plot.new()
for(i in 1:length(ls)){
  # from PAVO
  achrocol <- "grey"
  achrosize <- 1
  par(mar = c(1,1,2,2), bg = "white", pty = "s")
  vert <- data.frame(
    x = c(0, -1 / sqrt(2), 1 / sqrt(2)),
    y = c(sqrt(2) / sqrt(3), -sqrt(2) / (2 * (sqrt(3))), -sqrt(2) / (2 * (sqrt(3))))
  )
  
  arg <- list()
  arg$xlim <- c(-1 / sqrt(2), 1 / sqrt(2))
  arg$ylim <- c(-sqrt(2) / (2 * (sqrt(3))), sqrt(2) / sqrt(3))
  arg$x <- ls[[i]]$a$x
  arg$y <- ls[[i]]$a$y
  arg$xlab <- ""
  arg$ylab <- ""
  arg$bty <- "n"
  arg$axes <- FALSE
  
  do.call(plot, c(arg, type = "n"))
  polygon(vert, lwd = 1, lty = 1, border = "black")
  points(x = 0, y = 0, pch = 15, col = achrocol, cex = achrosize)
  
  # remove plot-specific args, add points after the stuff is drawn
  arg[c(
    "type", "xlim", "ylim", "log", "main", "sub", "xlab", "ylab",
    "ann", "axes", "frame.plot", "panel.first", "panel.last", "asp"
  )] <- NULL
  
  points(x = ls[[i]]$b$x, y = ls[[i]]$b$y,
         col = rgb(r = .1, g = 0.1, b = .1, alpha = 0.2), pch = 19)
  points(x = ls[[i]]$c$x, y = ls[[i]]$c$y,
         col = rgb(r = .3, g = 0.3, b = .3, alpha = 0.2), pch = 19)
  points(x = ls[[i]]$a$x, y = ls[[i]]$a$y,
         col = rgb(r = .8, g = 0.1, b = .1, alpha = 0.1), pch = 19)
  
  mtext(names(ls[i]), side = 3, line = -2, outer = TRUE)
  }
dev.off()
savep

points(x = ls[[1]]$a$x, y = ls[[1]]$a$y,
       col = "black", pch = 19)
points(x = ls[[1]]$b$x, y = ls[[1]]$b$y,
       col = rgb(r = 1, g = 0.2, b = 0, alpha = 0.1), pch = 19)
points(x = ls[[1]]$c$x, y = ls[[1]]$c$y,
       col = rgb(r = 0, g = 0.2, b = 1, alpha = 0.1), pch = 19)
points()
pch = 0:18

lrs <- unique(df$substrate) 
lrsList <- setNames(vector("list", length(lrs)), lrs)

lrsList <- lapply(lrs, function(x) df %>%
                 select(swMean, mwMean, lwMean, 
                        lumMean,pop_spp,substrate) %>%
                 filter(substrate %in% x))
lrsList <- setNames(lrsList, lrs)
tri.lrs <- lapply(lrsList, colspace, space = "tri")
par(mfrow = c(1,3))

lapply(names(tri.lrs), function(x) plot(tri.lrs[[x]], pch = 21, cex = 1, main = x))

