logit2 <- function(x) { log2(x) - log2(1-x) }

ilogit2 <- function(x) { 2^(x) / (1+2^(x)) }

.default.450k.annotation <- "ilmn.v1.2"
