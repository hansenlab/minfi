logit2 <- function(x) { log2(x) - log2(1-x) }

ilogit2 <- function(x) { exp(x) / (1+exp(x)) }

