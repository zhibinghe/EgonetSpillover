
#' @title Sample size (number of networks) calculation for a series of hypothesis tests of the Average Individual Effect (AIE) and
#' the Average Spillover Effect (ASpE)
#' @description These hypothesis tests include the direct tests for AIE (HIE) and ASpE (HSpE) respectively, the conjunctive (HISpJ) and dis-conjunctive (HISpC) test, and the test of overall effect (HOE).
#' @details
#' Let \eqn{\tau} and \eqn{\delta} be the AIE and ASpE respectively.
#'
#' \itemize{
#' \item{HIE:} The HIE test is used for the purpose with \eqn{H_0: \tau=0} against the alternative hypothesis \eqn{H_a: \tau \neq 0}.
#' A two-sided Z-test statistic \eqn{T_{\tau} = \sqrt{K} (\hat{\tau}/\hat{\sigma}_{\tau})}, where \eqn{\sigma^2_{\tau}} is the asymptotic variance of \eqn{\hat{\tau}}
#' (see equation (5) in the reference paper).
#' \item{HSpE:} The HSpE test is used for the purpose with
#' \eqn{H_0: \delta=0} against the alternative hypothesis \eqn{H_a: \delta \neq 0}.
#' A two-sided Z-test statistic \eqn{T_{\delta} = \sqrt{K} (\hat{\delta}/\hat{\sigma}_{\delta})}, where \eqn{\sigma^2_{\delta}} is the asymptotic variance of \eqn{\hat{\delta}}
#' (see equation (6) in the reference paper).
#' \item{HISpJ:} The HISpJ test is used for the purpose with
#' \eqn{H_0: \tau=0\  \text{and}\  \delta=0} against the alternative hypothesis \eqn{H_a: \tau \neq 0 \ \text{or}\ \delta \neq 0}.
#' A Wald test statistic \eqn{Q_J = K\hat{\bm{\theta}}^T \Sigma_{\tau\delta}^{-1} \hat{\bm{\theta}}}
#' is used and it follows a non-central \eqn{\chi^2} distribution (see equation (7) in the reference paper).
#' \item{HISpC:} The HISpC test is used for the purpose with
#' \eqn{H_0: \tau=0\  \text{or}\  \delta=0} against the alternative hypothesis \eqn{H_a: \tau \neq 0 \ \text{and}\ \delta \neq 0}.
#' A bivariate test statistic \eqn{Q_C = (T_{\tau},T_{\delta})} is used and it follows a bivariate normal distribution (see equation (8) in the reference paper).
#' \item{HOE:} The HOE test is used for the purpose with
#' \eqn{H_0: (\tau + n\delta)/(n+1)=0} against the alternative hypothesis \eqn{H_a: (\tau + n\delta)/(n+1) \neq 0}.
#' A two-sided Z-test statistic \eqn{T_{O} = \sqrt{K} \hat{\sigma}_{O} (\hat{\tau} + n\hat{\delta})/(n+1) }, where \eqn{\sigma^2_{O}} is the asymptotic variance of \eqn{(\hat{\tau} + n\hat{\delta})/(n+1)}
#' (see equation (9) in the reference paper).
#' }
#'
#' @param power size of power
#' @param tau Average Individual Effect (AIE)
#' @param delta Average Spillover Effect (ASpE)
#' @param rho_y Intra-class Correlation (ICC)
#' @param p random treatment allocation probability
#' @param n network size (social network members), in the basic ENRT design we assume the network size is equal for all networks
#' @param var_y variance of the outcome
#' @param alpha significant level
#' @param test choose which tests should be presented, by default, all the five tests (HIE, HSpE, HISpJ, HISpC and HOE) are presented
#' @importFrom mvtnorm pmvnorm
#' @return A vector with the following elements:
#' \itemize{
#' \item{K_HIE:} The minimal required number of networks for testing AIE
#' \item{K_HSpE:} The minimal required number of networks for testing ASpE
#' \item{K_HISpJ:} The minimal required number of networks for joint test of AIE and ASpE
#' \item{K_HISpC:} The minimal required number of networks for conjunctive test of AIE and ASpE
#' \item{K_HOE:} The minimal required number of networks for testing the overall effect
#' }
#' @export
#' @references  Design of egocentric network-based studies to estimate causal effects under interference (\url{https://arxiv.org/pdf/2308.00791})
#' @seealso \code{\link{n.HISp}}
#' @examples
#' alpha <- 0.05 # significant level
#' power <- 0.80 # power
#' tau <- 1  # effect size for testing the Average Individual Effect (AIE)
#' delta <- 1 # effect size for testing the Average Spillover Effect (ASpE)
#' n <- c(1, 2, 5,10) # network size
#' ICC <- seq(0.1, 0.9, by=0.1) # intra-class correlation
#' var_y <- 1 # variance of outcome
#' K.HISp(power=power, tau=tau, delta=delta, rho_y=0.1, p=0.5, n=5, var_y=var_y, alpha=alpha)
#' df <- expand.grid(ICC, n); colnames(df) <- c("ICC", "n")
#' out <- cbind(df, t(apply(df, 1, function(x) K.HISp(power=power, tau=tau, delta=delta, rho_y=x[1], p=0.5, n=x[2], var_y=var_y, alpha=alpha), simplify=TRUE)))
#' print(out)
K.HISp = function(power, tau, delta, rho_y, p, n, var_y, alpha, test=c('HIE','HSpE','HISpJ','HISpC','HOE')) {
  var_u <- rho_y*var_y # variance of random effect
  var_e <- var_y - var_u # variance of random noise
  var_tau <- var_y*(n*(1-p)*(1-rho_y )+(1+rho_y*n))/((n+1)*(1-p)*p) # variance of AIE
  var_delta <- var_y*((1-p)*(1-rho_y )+n*(1+n*rho_y ))/(n*(n+1)*(1-p)*p) # variance of ASpE
  var_td <- var_y*(p*(1+n*rho_y)+(1-p)*(n+1)*rho_y)/((n+1)*p*(1-p)) # covariance of AIE and ASpE
  #
  z_a <- qnorm(1-alpha/2, mean = 0, sd = 1) # quantile
  z_p <- qnorm(power, mean = 0, sd = 1)
  v <- qchisq(1-alpha, df=2, lower.tail = TRUE, log.p = FALSE)
  # KD: HIE
  # KS: HSpE
  # KDS: HISpJ
  # KDSC: HISpC
  # KO: overall effect test
  KD = KS = KDS = KDSC = KO <- NA
  if ('HIE' %in% test) KD <- var_tau * (z_a + z_p)^2/(tau^2)
  if ('HSpE' %in% test) KS <- var_delta * (z_a+z_p)^2/(delta^2)
  #
  v <- 9.633  # qchisq(0.95, 2, lower.tail = TRUE, log.p = FALSE) #? how to find the non-centrality parameter
  if ('HISpJ' %in% test) KDS <- (v*var_y*(1+n*rho_y))/((1-p)*p*(tau^2 + n*(delta^2)))
  #
  var_O <- var_y * (1 + n*rho_y)/((1+n)*(1-p)*p)
  if ('HOE' %in% test) KO <-  var_O * (z_a+z_p)^2/(tau+n*delta)^2 * (1+n)^2
  #
  f = function(K) {
    # power function for Conjunctive test, see eq (8) in the paper
    z_a2 <- qnorm(alpha/2, mean = 0, sd = 1)
    u_t <- sqrt(K)*tau/sqrt(var_tau)
    u_d <- sqrt(K)*delta/sqrt(var_delta)
    corr_td <- var_td/sqrt(var_tau*var_delta)
    Omega <- matrix(c(1, corr_td, corr_td, 1), 2, 2)
    intg1 <- pmvnorm(lower=c(z_a,z_a), upper=Inf, mean=c(u_t,u_d), Omega)
    intg2 <- pmvnorm(lower=c(z_a,-Inf), upper=c(Inf,z_a2), mean=c(u_t,u_d), Omega)
    intg3 <- pmvnorm(lower=c(-Inf,z_a), upper=c(z_a2,Inf), mean=c(u_t,u_d), Omega)
    intg4 <- pmvnorm(lower=c(-Inf,-Inf), upper=c(z_a2,z_a2), mean=c(u_t,u_d), Omega)
    return(intg1 + intg2 + intg3 + intg4 - power)
  }
  if ('HISpC' %in% test) KDSC <- suppressWarnings(rootSolve::multiroot(f, start=1, positive = TRUE, maxiter = 1000)$root)
  return(c(K_HIE=ceiling(KD), K_HSpE=ceiling(KS), K_HISpJ=ceiling(KDS), K_HISpC=ceiling(KDSC), K_HOE=ceiling(KO)))
}


#' @title Sample size (network size) calculation for a series of hypothesis tests of the Average Individual Effect (AIE) and
#' the Average Spillover Effect (ASpE)
#' @description These hypothesis tests include the direct tests for AIE (HIE) and ASpE (HSpE) respectively, the conjunctive (HISpJ) and dis-conjunctive (HISpC) test, and the test of overall effect (HOE).
#' @details
#' Let \eqn{\tau} and \eqn{\delta} be the AIE and ASpE respectively.
#'
#' \itemize{
#' \item{HIE:} The HIE test is used for the purpose with \eqn{H_0: \tau=0} against the alternative hypothesis \eqn{H_a: \tau \neq 0}.
#' A two-sided Z-test statistic \eqn{T_{\tau} = \sqrt{K} (\hat{\tau}/\hat{\sigma}_{\tau})}, where \eqn{\sigma^2_{\tau}} is the asymptotic variance of \eqn{\hat{\tau}}
#' (see equation (5) in the reference paper).
#' \item{HSpE:} The HSpE test is used for the purpose with
#' \eqn{H_0: \delta=0} against the alternative hypothesis \eqn{H_a: \delta \neq 0}.
#' A two-sided Z-test statistic \eqn{T_{\delta} = \sqrt{K} (\hat{\delta}/\hat{\sigma}_{\delta})}, where \eqn{\sigma^2_{\delta}} is the asymptotic variance of \eqn{\hat{\delta}}
#' (see equation (6) in the reference paper).
#' \item{HISpJ:} The HISpJ test is used for the purpose with
#' \eqn{H_0: \tau=0\  \text{and}\  \delta=0} against the alternative hypothesis \eqn{H_a: \tau \neq 0 \ \text{or}\ \delta \neq 0}.
#' A Wald test statistic \eqn{Q_J = K\hat{\bm{\theta}}^T \Sigma_{\tau\delta}^{-1} \hat{\bm{\theta}}}
#' is used and it follows a non-central \eqn{\chi^2} distribution (see equation (7) in the reference paper).
#' \item{HISpC:} The HISpC test is used for the purpose with
#' \eqn{H_0: \tau=0\  \text{or}\  \delta=0} against the alternative hypothesis \eqn{H_a: \tau \neq 0 \ \text{and}\ \delta \neq 0}.
#' A bivariate test statistic \eqn{Q_C = (T_{\tau},T_{\delta})} is used and it follows a bivariate normal distribution (see equation (8) in the reference paper).
#' \item{HOE:} The HOE test is used for the purpose with
#' \eqn{H_0: (\tau + n\delta)/(n+1)=0} against the alternative hypothesis \eqn{H_a: (\tau + n\delta)/(n+1) \neq 0}.
#' A two-sided Z-test statistic \eqn{T_{O} = \sqrt{K} \hat{\sigma}_{O} (\hat{\tau} + n\hat{\delta})/(n+1) }, where \eqn{\sigma^2_{O}} is the asymptotic variance of \eqn{(\hat{\tau} + n\hat{\delta})/(n+1)}
#' (see equation (9) in the reference paper).
#' }
#'
#' @param power size of power
#' @param tau Average Individual Effect (AIE)
#' @param delta Average Spillover Effect (ASpE)
#' @param rho_y Intra-class Correlation (ICC)
#' @param p random treatment allocation probability
#' @param K number of networks (number of index participants in the ENRT design)
#' @param var_y variance of the outcome
#' @param alpha significant level
#' @inheritParams K.HISp
#' @importFrom mvtnorm pmvnorm
#' @return A vector with the following elements:
#' \itemize{
#' \item{n_HIE:} The minimal required network size for testing AIE
#' \item{n_HSpE:} The minimal required network size for testing ASpE
#' \item{n_HISpJ:} The minimal required network size for joint test of AIE and ASpE
#' \item{n_HISpC:} The minimal required network size for conjunctive test of AIE and ASpE
#' \item{n_HOE:} The minimal required network size for testing the overall effect
#' }
#' @export
#' @references  Design of egocentric network-based studies to estimate causal effects under interference (\url{https://arxiv.org/pdf/2308.00791})
#' @seealso \code{\link{K.HISp}}
#' @examples
#' alpha <- 0.05 # significant level
#' power <- 0.80 # power
#' tau <- 1  # effect size for testing the Average Individual Effect (AIE)
#' delta <- 1 # effect size for testing the Average Spillover Effect (ASpE)
#' K <- c(20, 30) # network size
#' ICC <- seq(0.1, 0.8, by=0.1) # intra-class correlation
#' var_y <- 1 # variance of outcome
#' n.HISp(power=power, tau=tau, delta=delta, rho_y=0.1, p=0.5, K=30, var_y=var_y, alpha=alpha)
#' df <- expand.grid(ICC, K); colnames(df) <- c("ICC", "K")
#' out <- cbind(df, t(apply(df, 1, function(x) n.HISp(power=power, tau=tau, delta=delta, rho_y=x[1], p=0.5, K=x[2], var_y=var_y, alpha=alpha), simplify=TRUE)))
#' print(out)
n.HISp = function(power, tau, delta, rho_y, p, K, var_y, alpha, test=c('HIE','HSpE','HISpJ','HISpC','HOE')) {
  var_u <- rho_y*var_y # variance of random effect
  var_e <- var_y - var_u # variance of random noise
  z_a <- qnorm(1-alpha/2, mean = 0, sd = 1) # quantile
  z_p <- qnorm(power, mean = 0, sd = 1)
  v <- qchisq(1-alpha, df=2, lower.tail = TRUE, log.p = FALSE)
  # KD: HIE
  # KS: HSpE
  # KDS: HISpJ
  # KDSC: HISpC
  # KO: overall effect test
  n_HIE = n_HSpE = n_HISpJ = n_HISpJ = n_HOE <- NA
  KD.n = function(n) {
    return(K - var_y*(n*(1-p)*(1-rho_y)+(1+n*rho_y))*(z_a+z_p)^2/((n+1)*p*(1-p)*tau^2))
    #this is eq(5)
  }
  if ('HIE' %in% test) n_HIE <- suppressWarnings(rootSolve::multiroot(KD.n, start=1, positive = TRUE, maxiter = 1000)$root)
  #
  KS.n = function(n) {
    return(K - var_y*((1-p)*(1-rho_y)+n*(1+n*rho_y))*((z_a+z_p)^2)/(n*(n+1)*p*(1-p)*delta^2))
    # this is the eq(6) in the paper
  }
  if ('HSpE' %in% test) n_HSpE <- suppressWarnings(rootSolve::multiroot(KS.n, start=1, positive = TRUE, maxiter = 1000)$root)
  #
  v <- 9.633  # qchisq(0.95, 2, lower.tail = TRUE, log.p = FALSE) #? how to find the non-centrality parameter
  KDS.n = function(n) {
    # ?? confirm with Junhan about this magic number
    return(K - v*var_y*(1+n*rho_y)/(p*(1-p)*(tau^2+n*delta^2)))
    # this is eq(7) in the paper
  }
  if ('HISpJ' %in% test) n_HISpJ <- suppressWarnings(rootSolve::multiroot(KDS.n, start=1, positive = TRUE, maxiter = 1000)$root)
  #
  KO.n = function(n) {
    return(K - var_y*(1+n*rho_y)*((z_a+z_p)^2)*(n+1) / (p*(1-p)*((tau+n*delta)^2)) )
    # equation (9) in the paper
  }
  if ('HOE' %in% test) n_HOE <- suppressWarnings(rootSolve::multiroot(KO.n, start=1, positive = TRUE, maxiter = 1000)$root)
  #
  f = function(n){
    # this is the power function eq(8)
    z_a <- qnorm(1-alpha/2, mean = 0, sd = 1)
    z_a2 <- qnorm(alpha/2, mean = 0, sd = 1)
    var_td <- var_y*(p*(1+n*rho_y)+(1-p)*(n+1)*rho_y)/((n+1)*p*(1-p))
    var_tau <- var_y*(n*(1-p)*(1-rho_y)+(1+rho_y*n))/((n+1)*(1-p)*p)
    var_delta <- var_y*((1-p)*(1-rho_y)+n*(1+n*rho_y))/(n*(n+1)*(1-p)*p)
    u_t <- sqrt(K)*tau/sqrt(var_tau)
    u_d <- sqrt(K)*delta/sqrt(var_delta)
    corr_td <- var_td/sqrt(var_tau*var_delta)
    Omega <- matrix(c(1,corr_td,corr_td,1),2,2)
    intg1 <- pmvnorm(lower=c(z_a,z_a),upper=Inf,mean=c(u_t,u_d),Omega)
    intg2 <- pmvnorm(lower=c(z_a,-Inf),upper=c(Inf,z_a2),mean=c(u_t,u_d),Omega)
    intg3 <- pmvnorm(lower=c(-Inf,z_a),upper=c(z_a2,Inf),mean=c(u_t,u_d),Omega)
    intg4 <- pmvnorm(lower=c(-Inf,-Inf),upper=c(z_a2,z_a2),mean=c(u_t,u_d),Omega)
    return(intg1+intg2+intg3+intg4 - power)
  }
  if ('HISpC' %in% test) n_HISpC <- suppressWarnings(rootSolve::multiroot(f, start=1, positive = TRUE, maxiter = 1000)$root)
  return(c(n_HIE=ceiling(n_HIE), n_HSpE=ceiling(n_HSpE), n_HISpJ=ceiling(n_HISpJ), n_HISpC=ceiling(n_HISpC), n_HOE=ceiling(n_HOE)))
}

