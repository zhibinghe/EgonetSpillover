#' @title Noncentral Dunnett's test distribution
#' @description
#' For the purpose of checking, should be hidden in the final version
#'
#' @param q vector of quantiles \eqn{q \in R}
#' @param p vector of probabilities \eqn{(0, 1)}
#' @param rho vector of correlations, with length equal or greater than \eqn{1}
#' @param N vector size to be simulated, with \eqn{N > 1}
#' @param nu degrees of freedom \eqn{\nu > 0}
#' @param n number of points of the gaussian quadrature \eqn{n > 2}
#' @param delta vector of noncentrality parameter. Must be of the same length of \code{rho}
#' @param two.sided if true (default) the two-sided distribution is considered, else the one-sided distribution is considered
#' @details Assumes n = 32 as default value for dNCDun,  pNCDun and qNCDun. The \code{nu} can be a finite real number or an infinity real number. The default value of \code{nu} is \code{Inf} in \code{rNCDun}. When \code{nu=1}, the convergence of the routines requires \code{n>200} points in the Gaussian quadrature to obtain the desired result  accurately. The cumulative distribution  function of the noncentral unilateral Dunnett's test statistic with finite degrees of freedom \eqn{\nu}  is
#'     \deqn{F(q; r, \nu, \bm{\rho}, \bm{\delta})= \displaystyle \int_0^\infty \int_{-\infty}^\infty \phi(y) \prod_{j=1}^r \Phi\left(\frac{\sqrt{\rho_j} y +  x q-\delta_j}{\sqrt{1-\rho_j}}\right)  f(x;\nu)dy dx,}
#'    where \eqn{\bm{\rho}=[\rho_1, \rho_2, \ldots, \rho_r]^\top} is the correlation vector, \eqn{\bm{\delta}=[\delta_1, \delta_2, \ldots, \delta_r]^\top} is the vector of noncentrality parameter, \eqn{q} is the quantile of  the noncentral unilateral Dunnett's test distribution, \eqn{r} is the numbers of means (or sample size) and \eqn{\nu} is the degrees of freedom of a independent chi-square variable in the studentized process.The \eqn{f(x;\nu)} probability density function is given by
#'    \deqn{f(x; \nu)= \frac{ \nu^{\nu/2} }{\Gamma(\nu/2)2^{\nu/2-1}} x^{\nu-1} e^{-\nu x^2/2}, \quad x \ge 0.}
#'    The cumulative distribution  function of the noncentral unilateral Dunnett's test statistic with infinity degrees of freedom is
#'    \deqn{F(q; r, \nu=\infty, \bm{\rho}, \bm{\delta})=  \int_{-\infty}^\infty \phi(y) \prod_{j=1}^r \Phi\left(\frac{\sqrt{\rho_j} y + q-\delta_j }{\sqrt{1-\rho_j}}\right) dy.}
#'    The cumulative distribution  function of the noncentral bilateral  Dunnett's test statistic with finite degrees of freedom \eqn{\nu} is
#'    \deqn{F(q; r, \nu, \bm{\rho}, \bm{\delta})= \int_0^\infty \int_{-\infty}^\infty \phi(y) \prod_{j=1}^r \left[\Phi\left(\frac{\sqrt{\rho_j} y +  x q-\delta_j}{\sqrt{1-\rho_j}}\right) - \Phi\left(\frac{\sqrt{\rho_j} y -  x q-\delta_j}{\sqrt{1-\rho_j}}\right)\right]  f(x;\nu)dy dx.}
#'    Finally, the cumulative distribution  function of the noncentral bilateral  Dunnett's test statistic with infinity degrees of freedom is
#'    \deqn{F(q; r, \nu=\infty, \bm{\rho}, \bm{\delta})= \int_{-\infty}^\infty \phi(y) \prod_{j=1}^r \left[\Phi\left(\frac{\sqrt{\rho_j} y + q -\delta_j}{\sqrt{1-\rho_j}}\right) - \Phi\left(\frac{\sqrt{\rho_j} y - q-\delta_j }{\sqrt{1-\rho_j}}\right)\right] dy.}
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{dNCDun} gives the density
#' \item{pNCDun} gives the cumulative distribution function
#' \item{qNCDun} gives the quantile function
#' }
#' @export
#' @examples
#' n <- 50
#' nu <- 9
#' rho <- c(0.5,0.5,0.5)
#' delta <- c(0,0,0)
#' q <- 2.30
#' p <- 0.95
#' pNCDun(q, nu, rho, delta, 32, TRUE)
#' dNCDun(q, nu, rho, delta, 32, TRUE)
#' qNCDun(p, nu, rho, delta, 16, TRUE)
#' q <- c(2.34, 4.50, 3.40)
#' p <- c(0.85, 0.95, 0.975)
#' nu <- c(Inf, 19, 15)
#' # unbalanced example
#' rho   <- c(0.23, 0.25, 0.27)
#' delta <- c(0, 0, 0) # central case
#' pNCDun(q, nu, rho, delta, 32, TRUE)
#' dNCDun(q, nu, rho, delta, 32, TRUE)
#' qNCDun(p, nu, rho, delta, 16, TRUE)
#'

dNCDun <- function(q, nu, rho, delta, n = 32, two.sided = TRUE)
{
  x <- GaussLegendre(n)
  nn <- length(q)
  if (nn == length(nu))  xx <- cbind(q, nu) else
    if (length(nu) == 1) xx <- cbind(q, rep(nu, times = nn))
  dtched <- function(xx)
  {
    if (two.sided==TRUE)
    {
      if (xx[2]==Inf) return(dNDBD(xx[1], rho, delta, n, x)) else
        return(dNDBDF(xx[1], rho, xx[2], delta, n, x))
    } else
    {
      if (xx[2]==Inf) return(dNDUD(xx[1], rho, delta, n, x)) else
        return(dNDUDF(xx[1], rho, xx[2], delta, n, x))
    }
  }
  d <- apply(xx, 1, dtched)
  return(d)
}

# cdf of Dunnett distribution
pNCDun <- function(q, nu, rho, delta, n = 32, two.sided = TRUE)
{
  x <- GaussLegendre(n)
  nn <- length(q)
  if (nn == length(nu))  xx <- cbind(q, nu) else
    if (length(nu) == 1) xx <- cbind(q, rep(nu, times = nn))
  dtched <- function(xx)
  {
    if (two.sided==TRUE)
    {
      if (xx[2]=="Inf") return(pNDBD(xx[1], rho, delta, n, x)) else
        return(pNDBDF(xx[1], rho, xx[2], delta, n, x))
    } else
    {
      if (xx[2]=="Inf") return(pNDUD(xx[1], rho, delta, n, x)) else
        return(pNDUDF(xx[1], rho, xx[2], delta, n, x))
    }
  }
  p <- apply(xx, 1, dtched)
  return(p)
}

# quantile of Dunnett distribution
qNCDun <- function(p, nu, rho, delta, n = 32, two.sided = TRUE)
{
  x <- GaussLegendre(n)
  nn <- length(p)
  if (nn == length(nu))  xx <- cbind(p, nu) else
    if (length(nu) == 1) xx <- cbind(p, rep(nu, times = nn))
  dtched <- function(xx)
  {
    if (two.sided==TRUE)
    {
      if (xx[2]==Inf) return(qNDBD(xx[1], rho, delta, n, x)) else
        return(qNDBDF(xx[1], rho, xx[2], delta, n, x))
    } else
    {
      if (xx[2]==Inf) return(qNDUD(xx[1], rho, delta, n, x)) else
        return(qNDUDF(xx[1], rho, xx[2], delta, n, x))
    }
  }
  d <- apply(xx, 1, dtched)
  return(d)
}


