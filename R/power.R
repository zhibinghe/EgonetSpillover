#' @title MCB standard power calculation based on given sample size and number of groups
#'
#' @param n the size of each group (currently we assume it is the same across all groups)
#' @param alpha significant level
#' @param k number of groups
#' @param w effect size in alternative hypothesis
#' @references Hsu, J. (1996). Multiple comparisons: theory and methods. CRC Press.
#' @return mcb power
#' @export
#' @examples
#' n.sam = 2:12
#' xx = sapply(n.sam, power.size.mcb, k=8, w = 5/3, alpha = 0.05)
#' plot(n.sam,xx,type="l",col="red", ylab="Power", xlab="sample size")
#' abline(h = 1-0.05, lty=2,col="grey")
#' abline(h=0.90,lty=2,col="blue")
#' # simulation example in Hsu's textbook on page-241
#' alpha = 0.10
#' n.sam = seq(2, 24, by=1)
#' ratio = seq(0.8,1.2,by=0.1)
#' dat.pwr = matrix(NA,nrow = length(n.sam), ncol = length(ratio))
#' for (j in 1:length(ratio)) dat.pwr[,j] = sapply(n.sam, power.size.mcb, k=5, w = ratio[j], alpha = alpha)
#' plot(lowess(n.sam,dat.pwr[,1],f = 0.1), type="l",col=1,
#'      xlab = "sample size", ylab = "Power")
#' abline(h = 1-alpha, lty = 2, col="grey")
#' for (j in 2:length(ratio)) {
#'   lines(lowess(n.sam,dat.pwr[,j],f = 0.1),col=j)
#' }
#' legend("bottomright",legend = paste0("ratio = ", ratio), col = 1:length(ratio),
#'        lty = rep(1,length(ratio)), bty="n")

power.size.mcb = function(n, alpha = 0.05, k, w){
  # rho: k-1 length of vector
  nu = k * (n-1)
  power.f = function(z,s,q,df) {
    (pnorm(z + sqrt(2) * q * s))^(k-1) * dnorm(z) * chi::dchi(s*sqrt(df),df) * sqrt(df)
  }
  d = nCDunnett::qNCDun(p = 1-alpha, nu = nu, rho = rep(0.5,k-1), delta = 0, two.sided = FALSE)
  u = w*sqrt(n/2)/d
  pracma::integral2(power.f, xmin=-6,xmax=6, ymin=0, ymax=u,q=d,df = nu)$Q
}



