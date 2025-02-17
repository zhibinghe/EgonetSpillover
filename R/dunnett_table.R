# Critical value calculation for multiple comparison with control (MCC) and multiple comparison with the best (MCB)
# q probability
# g number of groups
# tail if true (default) the two-sided distribution is considered, else the one-sided distribution is considered, see qmvt
# references Hsu, J. (1996). Multiple comparisons: theory and methods. CRC Press.
# return the critical value \eqn{c_{\alpha}}
# import mvtnorm
# examples
# library(mvtnorm) ## for multivariate t-distribution
# g.sizes <- c(2:9)
# critval.two <- critval.one <- numeric(length(g.sizes))
#
# for(j in seq_along(g.sizes)){
#   critval.two[j] <- qDunnett(q = 0.95, g.sizes[j], df = 10, tail = "both.tails")
#   critval.one[j] <- qDunnett(q = 0.99, g.sizes[j], df = 10, tail = "lower.tail")
# }
# round(critval.two, 2) # compare with Hsu's textbook, page 250
# round(critval.one, 2) # compare with Hsu's textbook, page 254
# # Use qNCDun function
# quantiles <- sapply(g.sizes, function(x) qNCDun(p = 0.95, nu = 10, rho = rep(0.5, x - 1), delta = 0))
# round(quantiles, 2)
qDunnett <- function(q, g, df, tail){

  corr <- matrix(0.5, nrow = g - 1, ncol = g - 1) ## implicitly assumes *balanced* design!
  diag(corr) <- 1

  if(length(g) > 1)
    stop("only scalar g allowed")
  if(length(q) > 1)
    stop("only scalar q allowed")

  qmvt(q, tail = tail, df = df, corr = corr)$quantile
}


