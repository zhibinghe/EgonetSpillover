################################################################################
############# Ego-network design with multivariate X
############# Sample size calculation for spillover effect of heterogeneity
################################################################################

# Compute Variance matrix from correlation matrix

cor2cov = function(R, s) diag(s) %*% R %*% diag(s)

# Compute Variance matrix of categorical random variable
# g: a vector of probabilities for each category, the sum is 1
# return variance matrix

cov.cat = function(g) -outer(g, g) + diag(g)

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

#' Sample size calculation for Wald test with continuous covariates
#'
#' @param p Bernoulli probability for treatment of index participants
#' @param rho ICC correlation
#' @param n network size (exclude ego node), suppose all networks have the same number of members
#' @param m number of covariates
#' @param sigma2_y variance of outcome variable
#' @param sigma2_x variance matrix of covariates
#' @param Delta a vector of values in alternative hypothesis
#' @param alpha Type I error, 0.05 by default
#' @param lambda Type II error
#' @param Kmax maximum iteration times
#' @param message logical value indicates whether running procedure is reported
#'
#' @return minimal required sample size
#' @import dplyr
#' @export
#' @examples
#' library(ggplot2)
#' m <- 3
#' R <- matrix(0.0, m, m); diag(R) <- 1
#' s <- rep(1, m)
#' Delta <- rep(0.2, m) # repeat case
#' # Delta <- c(0.2, rep(0, m-1)) #unrepeat case
#' sam_size.Wald.cnt(p = 0.3, rho=0.1, m=m, sigma2_x=cor2cov(R,s), n=2, Delta=Delta, message=TRUE) # simple demo

#' ## simulation for various parameters combination
#' var_n <- c(2,3,5,10)
#' #var_p <- c(0.3,0.5,0.7,0.9)
#' var_p <- 0.5
#' var_rho <- seq(0.1, 0.9, 0.1)
#' outK <- list(length = length(var_n))

#' for(i in 1:length(var_n)) {
#'   svMisc::progress(i, length(var_n))
#'   temp <- data.frame(matrix(NA, nrow=length(var_rho), ncol=length(var_p) + 1))
#'   temp[,1] <- var_rho
#'   for(j in 1:length(var_p)) {
#'     temp[,j+1] <- sapply(var_rho, sam_size.Wald.cnt, p=var_p[j], n=var_n[i], m=m, sigma2_x=cor2cov(R,s), Delta=Delta) #
#'   }
#'   outK[[i]] <- temp
#' }

#' outK %>% dplyr::bind_rows(.id = "id") %>% tidyr::pivot_longer(-c(id,X1)) %>%
#'   ggplot(aes(x = X1, y = value, color = name, group = name)) +
#'   geom_smooth(method="loess", span=0.7, se=F) +
#'   facet_wrap(~id, labeller=labeller(id = c("1" = "n = 2",
#'                                            "2" = "n = 3",
#'                                            "3" = "n = 5",
#'                                            "4" = "n = 10"))) +
#'   labs(y = paste0("Required number of networks"), x = expression(rho[Y])) +
#'   scale_color_discrete(name="covariates", labels=c(3,4))
#'
sam_size.Wald.cnt = function(p, rho, n, m=1, sigma2_y=1, sigma2_x=NULL, Delta,
                             alpha=0.05, lambda=0.2, Kmax=1000, message=FALSE) {
  # sigma2_x is identity matrix if missing
  if (is.null(sigma2_x)) {
    if (m == 1) sigma2_x <- 1
    if (m > 1) sigma2_x <- diag(rep(1, m))
  }
  # 1D case -- closed-form solution
  if (m == 1) {
    # variance of estimator
    sigma2_deltax <- sigma2_y/sigma2_x * (n+1) * (1-rho) * (1+n*rho) /(n * p * ((n+1) - n*p*(1-rho)))
    sampsize <- sigma2_deltax/Delta^2 * (qnorm(1-alpha/2) + qnorm(1-lambda))^2
  }
  # multiple dimensional case -- Wald test
  if (m > 1) {
    if (rho == 1) k <- 0 # save time
    else{
      # variance matrix of estimators
      sigma2_deltax <- sigma2_y * (n+1) * (1-rho) * (1+n*rho) /(n*p * ((n+1)-n*p*(1-rho))) * solve(sigma2_x)
      chisqf = function(x, KK) {
        ncp <- KK * t(Delta) %*% solve(sigma2_deltax) %*% Delta
        return(dchisq(x, df=m, ncp=ncp))
      }
      # starting values
      k <- 1 ; intg <- 0
      while (intg < 1-lambda & k <= Kmax) {
        if(message) { if (k %% 10 == 0) cat(paste0("now processing K = ", k),"\n") }
        intg <- integrate(chisqf, lower=qchisq(1-alpha, m), upper=Inf, KK=k)$value
        k <- k + 1
      }
    }
    sampsize <- k
  }
  return(sampsize)
}

#' Sample size calculation for Wald test with categrical covariates
#'
#' @param p Bernoulli probability for treatment of index participants
#' @param rho ICC correlation
#' @param n network size (exclude ego node), suppose all networks have the same number of members
#' @param m number of dummy variables, number of categories - 1
#' @param g a vector of probabilities for each category, the sum is 1
#' @param sigma2_y variance of outcome variable
#' @param sigma2_x variance matrix of covariates
#' @param Delta a vector of values in alternative hypothesis
#' @param alpha Type I error, 0.05 by default
#' @param lambda Type II error
#' @param Kmax maximum iteration times
#' @param message logical value indicates whether running procedure is reported
#'
#' @return minimal required sample size
#' @export
#' @import dplyr
#' @examples
#' # example code
#' m <- 5 # 3,4,5
#' Delta <- rep(0.5, m)
#' # Delta <- c(0.5, rep(0, m-1))
#' g <- rep(1/(m+1), m+1)
#' sam_size.Wald.cat(p=0.3, rho=0.1, n=5, g=g, m=m, sigma2_x=cov.cat(g), Delta=Delta, message=TRUE)

#' var_n <- c(2,3,5,10)
#' var_p <- 0.5
#' var_rho <- seq(0.1, 0.9, 0.1)
#' outH <- list(length = length(var_n))
#' for(i in 1:length(var_n)) {
#'   svMisc::progress(i, length(var_n))
#'   temp <- data.frame(matrix(NA, nrow = length(var_rho), ncol=length(var_p) + 1))
#'   temp[,1] <- var_rho
#'   for(j in 1:length(var_p)) {
#'     temp[,j+1] = sapply(var_rho, sam_size.Wald.cat, g=g, p=var_p[j], n=var_n[i], sigma2_x=cov.cat(g), Delta=Delta)
#'   }
#'   outH[[i]] <- temp
#' }
#' outH %>% dplyr::bind_rows(.id = "id") %>% tidyr::pivot_longer(-c(id,X1)) %>%
#'   ggplot(aes(x = X1, y = value, color = name, group = name)) + # y = value * 8
#'   geom_smooth(method="loess", span=0.9, se=F) +
#'   labs(y = paste0("Required number of networks"), x = expression(rho[Y])) +
#'   scale_color_discrete(name="categories", labels=c("H=3","H=4","H=5","H=6")) +
#'   scale_linetype_manual(values = c("X2"=1,"X3"=1,"X4"=2,"X5"=2)) +
#'   #scale_linetype_manual(name = "categories", values = c(2,2,1,1))
#'   facet_wrap(~id, labeller = labeller(id = c("1" = "n = 2",
#'                                              "2" = "n = 3",
#'                                              "3" = "n = 5",
#'                                              "4" = "n = 10")))
#'
sam_size.Wald.cat = function(p, rho, n, g, sigma2_y=1, sigma2_x=NULL, Delta=0.5,
                             alpha=0.05, lambda=0.2, Kmax=1000, message=FALSE, maineff=FALSE) {
  # Variance matrix of estimators
  # return variance matrix of for X(-i) # i is the reference level
  m <- length(g) - 1 # number of dummy variables, number of categories -1
  Vi <- function(n_net, g) {
    if (maineff) temp <- cov.delta_maineff(n_net, n_memb=n, rho, p, g, sigma2_y)
    if (!maineff) temp <- cov.delta(n_net, n_memb=n, rho, p, g, sigma2_y)
    Ci = cbind(diag(rep(1, m)), rep(1, m))
    return(Ci %*% temp %*% t(Ci))
  }
  if (rho == 1) k <- 0
  else {
    chisqf = function(x, KK) {
      ncp <- KK * t(Delta) %*% solve(Vi(KK,g)) %*% Delta
      dchisq(x, df=m, ncp = ncp)
    }
    # starting value
    k <- 1 ; intg <- 0
    while (intg < 1-lambda & k <= Kmax) {
      if(message) {if (k %% 10 == 0) cat(paste0("now processing K = ", k),"\n")}
      intg <- integrate(chisqf, lower=qchisq(1-alpha, m), upper=Inf, KK=k)$value
      k <- k + 1
    }
  }
  sampsize <- k
  return(sampsize)
}


