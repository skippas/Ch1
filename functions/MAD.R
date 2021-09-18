# calculates Mean Absolute Deviation of a vector. Use inside of summarize
MAD <- function(x){sum(abs(x - mean(x)))*(1/length(x))}

