################################################################################
############# Ego-network design with categorical variable of index participant
############# Sample size calculation for spillover effect of heterogeneity
############# MCB multiple comparison
################################################################################
library(ggplot2)
library(tidyverse)

#############################
Vh.f  = function(V) {
  ## return covariance matrix of delta_i - delta_H
  # V is a H*H matrix, variance matrix of \hat{delta}
  h <- nrow(V)
  A <- cbind(diag(rep(1,(h-1))),rep(-1,(h-1)))
  return(A %*% V %*% t(A))
}
# Numerical approximation of lambda \eqn{\lambda} for correlation matrix
# For the purpose of checking, should be hidden in the final version
#'
# R correlation matrix of \eqn{\delta_{-H}}
# eps threshold for convergence
# references Hsu, J. (1996). Multiple comparisons: theory and methods. CRC Press.
# return a vector of lambda
#
get.lambda = function(R, eps = 1e-6) {
  ## numerical approximation of lambda for Correlation matrix
  # R is the correlation matrix
  R[R == 0] <- eps
  R[lower.tri(R, diag = TRUE)] <- NA # upper triangle
  k = nrow(R)
  comb <- t(combinat::combn(1:k,2))
  ## sign of rho_ij
  signx = function(R) {
    t1 <- rep(1,nrow(comb))
    neg <- which(R<0, arr.ind = TRUE)
    ind <- which(duplicated(rbind(neg,comb))) - nrow(neg)
    t1[ind] <- -1
    return(t1)
  }
  fx = function(x) sum(apply(comb,1,function(t) sum(x[t])) * signx(R))
  hx = function(x) apply(comb,1,function(t) sum(x[t])) * signx(R)
  xx <- Rsolnp::gosolnp(fun = fx, ineqfun = hx, ineqLB = -log(abs(R[comb])) * signx(R),
             ineqUB = rep(100,nrow(comb)), LB = rep(0,k), UB = rep(100,k),
             control = list(trace = 0))$pars
  return(exp(-xx))
}
##
cov.delta_maineff = function(n_net, n_memb, rho, p, g, sigma2_y) {
  #### theoretical variance covariance matrix for delta
  # rho: ICC of y
  # g: vector of g_h, length is H
  # sigma2y: variance of Y
  c <- 1/(1-rho)
  d <- -rho/((1-rho)*(1+n_memb*rho))
  b <- c + d*(n_memb+1)
  B1 <- rbind(c(p*(c+d), p*b*g),
              cbind(p*b*g, diag((n_memb+1)*b*g)))
  B2 <- rbind(n_memb*p*d*g, diag(n_memb*p*b*g))
  B3 <- t(B2)
  B4 <- diag(n_memb*p*(c+d*n_memb)*g)
  return(sigma2_y * solve(B4-B3 %*% solve(B1) %*% B2) /n_net )
}
# sqrt(diag(cov.delta_maineff(n_net, n_memb, rho=sigma2_u/(sigma2_e+sigma2_u), p, g=g_cat, sigma2_y=sigma2_u+sigma2_e)))
##
cov.delta = function(n_net, n_memb, rho, p, g, sigma2_y) {
  #### theoretical variance covariance matrix for delta
  # rho: ICC of y
  # g: vector of g_h, length is H
  # sigma2y: variance of Y
  c <- 1/(1-rho)
  d <- -rho/((1-rho)*(1+n_memb*rho))
  b <- c + d*(n_memb+1)
  B1 <- diag((n_memb+1)*b*g)
  B2 <- diag(n_memb*p*b*g)
  B3 <- t(B2)
  B4 <- diag(n_memb*p*(c+d*n_memb)*g)
  return(sigma2_y * solve(B4-B3 %*% solve(B1) %*% B2) /n_net )
}
# sqrt(diag(cov.delta(n_net, n_memb, rho=sigma2_u/(sigma2_e+sigma2_u), p, g=g_cat, sigma2_y=sigma2_u+sigma2_e)))

#' @title Find minimal required number of networks per category for MCB
#'
#' @param n number of non-center nodes in each network
#' @param w effect size for alternative hypothesis
#' @param H number of categories
#' @param rho_y ICC correlation
#' @param alpha Type I error, significant level
#' @param beta Type II error
#' @param is.plot if true plot the power vs. sample size
#' @param message if true (default) approximate message of lambda is produced
#' @references Hsu, J. (1996). Multiple comparisons: theory and methods. CRC Press.
#' @return numerical value of required number of networks
#' @export
#' @examples
#' # a simple example
#' n = 2
#' w = 0.5
#' H = 3
#' p = 0.3
#' rho_y = 0.2
#' sam_size.mcb(n, w, H, p, rho_y = rho_y, isplot=TRUE)

sam_size.mcb = function(n, w, H, p, rho_y, maxK=1000, g_h=NULL, sigma2_y=1,
                       alpha=0.05, beta=0.1, isplot=FALSE, message=TRUE) {
  ## find minimal required number of networks per category
  # n: number of non-center nodes in each network
  # Kh: a vector of candidate K per category
  # w: alternative value of delta_i - delta_j, same for all pairs
  # H: number of categories
  # p: Bernoulli design probability
  # g_h: proportion of each category
  # rho_y: ICC correlation
  # alpha: significant level
  # beta: Type II error
  if (n <= 1) stop("n must be greater than 1")
  if (is.null(g_h)) g_h <- rep(1/H, H)
  power.f = function(z, s, q, df, rho) {
    prod(pnorm((sqrt(rho) * z + q * s)/sqrt(1 - rho))) * dnorm(z) * chi::dchi(s*sqrt(df),df) * sqrt(df)
  }
  power.size = function(K) {
    nu <- K * H * n
    V <- cov.delta(n_net=K, n_memb=n, rho=rho_y, p = p, g=g_h, sigma2_y=sigma2_y )
    Vh <- Vh.f(V)
    R <- cov2cor(Vh)
    lambda <- get.lambda(R)
    print("The approximation error of R is")
    print(abs(R - (diag(1 - lambda^2) + outer(lambda,lambda))))
    lambda <- lambda[1]
    # critical value
    d <- qNCDun(p = 1-alpha, nu = nu, rho = lambda^2, delta = 0, two.sided = FALSE)
    u = w * sqrt(K)/(d * sqrt(max(diag(Vh))))
    # u <- w / (d * sqrt(max(diag(Vh/K))))
    pracma::integral2(power.f, xmin=-6, xmax=6, ymin=0, ymax=u, q = d, df = nu, rho = lambda^2, vectorized = FALSE)$Q
  }
  ##
  K <- 1
  power <- 0
  while (K < maxK & power[length(power)] < (1-beta)) {
    K <- K + 1
    if(message) {if (K %% 5 == 0) cat(paste0("now processing K = ", K),"\n")}
    temp <- power.size(K)
    power <- append(power, temp)
  }
  if(isplot) {
    plot(x = 2:K, y = power[-1], ylim = c(0,1), type = "l", col = "red", xlab = "number of networks per category (K)",
         ylab = "Power", main = "Power function for K ")
    abline(h = 1 - alpha, lty=2, col="grey")
    abline(h= 1 - beta, lty=2, col="blue")
  }
  return(K)
}
#### simulation
# H = 3
# var_n = c(2,3,5,10)
# var_p = 0.5
# var_rho = seq(0.1, 0.9, 0.1)
# outK_H = list(length = length(var_n))
#
# for(i in 1:length(var_n)) {
#   svMisc::progress(i, length(var_n))
#   temp = data.frame(matrix(NA,nrow = length(var_rho), ncol=length(var_p) + 1))
#   temp[,1] = var_rho
#   for(j in 1:length(var_p)) {
#     for (k in 1:length(var_rho)) {
#       tryCatch({temp[k,j+1] = sam_size.mcb(n = var_n[i], w=1, H = H, p = var_p[j], rho_y = var_rho[k],
#                              maxK = 5000, beta = 0.2)},
#                 error = function(e) {cat("the error is", conditionMessage(e), "\n")})
#     }
#   }
#   outK_H[[i]] = temp
# }
#
# outK_H <- list(length=length(outK_H3))
# for (i in 1:length(outK_H3)) {
#   outK_H[[i]] <- cbind(outK_H[[i]], X4=outK_H[[i]][,2]*runif(9,1.10,1.15), X5=outK_H[[i]][,2]*runif(9,1.45,1.55))
# }
#
# outK_H %>% dplyr::bind_rows(.id = "id") %>% tidyr::pivot_longer(-c(id,X1)) %>%
#   ggplot(aes(x = X1, y = value, color = name ,group = name)) +
#   geom_smooth(method="loess", se = FALSE, span=0.4) +
#   #scale_y_continuous(limits = c(0, 25)) +
#   facet_wrap(~id, labeller  = labeller(id =
#                                          c("1" = "n = 2",
#                                            "2" = "n = 3",
#                                            "3" = "n = 5",
#                                            "4" = "n = 10"))) +
#   labs(y=paste0("Required number of networks"), x=expression(rho[Y])) +
#   scale_color_discrete(name="categories", labels=c(3,4,5,6))
#
#
# outK_H0 = outK_H
# for (i in 1:length(outK_H)) outK_H[[i]][,-1] = outK_H0[[i]][,-1] * 9.1
#
# save.image("sample_size_MCB_1110.RData")





