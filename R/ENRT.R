
# add the network graph -- when network is small
#' @title Convert to the standard ego-network data structure based on users' data
#'
#' @param nodeid a vector of unique id for all network nodes
#' @param R a vector of logical values indicating if a node is an index participant, has values of either 0 or 1, must have the same length as \code{nodeid}
#' @param idx_net a vector showing which network each node belongs to, must have the same length as \code{nodeid}. Note that each network member can only be connected to one index participant
#' @param Z a vector of treatment status for all nodes, must have the length as \code{nodeid}
#' @param X a vector or a matrix of continuous or categorical covariates for each node
#' @param y a vector of outcome for each node, must have the same length as \code{nodeid}
#' @details
#' The standard ego-network dataframe should be organized and grouped by the network, and include the following required columns information:
#' \itemize{
#' \item{netid:} network id, such as \eqn{k=1,\dots,K}
#' \item{nodeid:} nodes in each network
#' \item{R:} indicator of index participant, showing if a node is the index participant of a network, ENRT assumes each network has and only has one index participant
#' \item{Z:} treatment status, indicating if the node is treated. ENRT assumes that only the index participant can be treated
#' \item{G:} number of treated social network neighbors
#' \item{X:} covariates, if there are multiple categorical variables, they must be integrated to one categorical variable.
#' \item{y:} outcome variable
#' }
#' @import dplyr
#' @return a standard dataframe of ego-network
#' @seealso \code{gen_egonet}
#' @examples
#' # random generate ENRT data
#' n_node <- 100
#' K <- 20
#' nodeid <- 1:n_node
#' netid <- 1:K
#' idx_net <- sample(netid, n_node, replace=TRUE)
#' # might give error, data does not satisfy the ENRT design
#' R <- sample(c(0,1), n_node, replace=TRUE)
#' Z <- sample(c(0,1), n_node, replace=TRUE)
#' df.ENRT(nodeid, R, idx_net, Z)
#' # ENRT design
#' df <- data.frame(nodeid=nodeid, idx_net=idx_net) %>% arrange(idx_net, by_group=TRUE)
#' # add R
#' df <- df %>% group_by(idx_net) %>% mutate(R= ifelse(row_number()==1, 1, 0))
#' # add Z
#' df <- df %>% mutate(Z= ifelse(R==1, sample(c(0,1), 1, prob=c(0.6,0.4), replace=TRUE), 0))
#' df.ENRT(df$nodeid, df$R, df$idx_net, df$Z, X=sample(LETTERS[1:4], n_node, replace=TRUE, prob=c(1/4,1/4,1/4,1/4)))

df.ENRT = function(nodeid, R, idx_net, Z, X=NA, y=NA) {
  df <- data.frame(nodeid=nodeid, R=R, netid=idx_net, Z=Z) %>% arrange(netid, desc(R), by_group=TRUE)
  # check if each network has only one index participant and only the index participants can be treated
  netRsum <- df %>% group_by(netid) %>% summarise(netsum = sum(R))
  if (any(netRsum$netsum != 1))
    stop(paste('These network have problems:', paste(netRsum$netid[netRsum$netsum != 1], collapse=" "), ', each network has and only has one index participant'))
  # check if only the index participants can be treated
  netRZ <- df %>% filter(R==0, Z==1)
  if (nrow(netRZ) !=0)
    stop(paste('These network have problems:', paste(netRZ$netid, collapse=" "), ', only the index participants can be treated'))
  # compute G
  # network id
  df$G <- 0 # initial
  n_netid <- df %>% group_by(netid) %>% summarise(n_net=n(), G=sum(Z))
  for (id in n_netid$netid) {
    df[df$netid == id, 'G'] <- unlist(c(0, rep(n_netid[n_netid$netid==id, 'G'], n_netid[n_netid$netid==id, 'n_net']-1)))
  }
  # add covariates and outcome
  df <- df %>% mutate(X=X, y=y)
  return(df)
}


#' @title Model spillover effect of ENRT data using mixed effect model
#' @details
#'  The following mixed effect model is applied to estimate the spillover effect of ENRT data
#' \deqn{Y_{ik}  = \sum_{h = 1}^{H} \zeta_{h} S_{kh} + \sum_{h = 1}^{H} \delta_{h} G_{ik}S_{kh} + u_k + \epsilon_{ik}, i = 2,\dots, n+1, k= 1,\dots, K,}
#' where \eqn{Y_{ik}} is the outcome for node \eqn{ik} in network \eqn{k}, \eqn{S_{kh}} is the category indicator for index participant in network \eqn{k} belonging to category \eqn{h},
#' \eqn{G_{ik}} is the number of treated network members for node \eqn{ik},
#' \eqn{\delta_h} represents the spillover effect of a network member whose index participant belongs to category \eqn{h}.
#' \eqn{u_k} is the random effect for networks and \eqn{\epsilon_{ik}} is the random noise.
#' @param data The standard ENRT data with categorical covariate X and outcome y, created from the function \code{df.ENRT}
#' @seealso \code{\link{df.ENRT}}, \code{\link[lmer]{lmer}}
#' @return # ...... the RMLE estimate ...
#' @export
#' @examples
#' # simulation
#' n_node <- 1000
#' K <- 100
#' nodeid <- 1:n_node
#' netid <- 1:K
#' # ENRT design
#' df <- data.frame(nodeid=nodeid, idx_net=idx_net) %>% arrange(idx_net, by_group=TRUE)
#' # add R
#' df <- df %>% group_by(idx_net) %>% mutate(R= ifelse(row_number()==1, 1, 0))
#' # add Z
#' df <- df %>% mutate(Z= ifelse(R==1, sample(c(0,1), 1, prob=c(0.6,0.4), replace=TRUE), 0)) # random treatment allocation
#' # add X
#' df$X <- sample(LETTERS[1:4], n_node, replace=TRUE, prob=c(0.1,0.2,0.3,0.4))
#' df.ENRT(df$nodeid, df$R, df$idx_net, df$Z, X=df$X)
#' # add y
#' data <- df.ENRT(df$nodeid, df$R, df$idx_net, df$Z, X=df$X)
#' S <- as.matrix(fastDummies::dummy_cols(df$X, remove_selected_columns=TRUE))
#' cat_neigh <- data$G * S
#' # parameter setup
#' zeta <- rep(1, 4)
#' delta <- rep(2, 4)
#' sigma2_u <- 36
#' sigma2_e <- 4
#' y <- as.vector(S %*% zeta) + as.vector(cat_neigh %*% delta) + rnorm(K, mean=0, sd=sqrt(sigma2_u)) +
#'  rnorm(n_node, mean=0, sd=sqrt(sigma2_e))
#' data_full <- df.ENRT(df$nodeid, df$R, df$idx_net, df$Z, X=df$X, y=y)
#' model <- est.ENRT.sp(data_full) # lmer model
#' summary(model)
#' vcov(model)[5:8,5:8]
#' # need to add Junhan's model for estimation of AIE and ASpE

est.ENRT.sp = function(data) {
  #
  Skh <- as.matrix(fastDummies::dummy_cols(data$X, remove_selected_columns=TRUE))
  inter_cat_neigh <- data$G * Skh
  model <- lme4::lmer(y ~ 0 + Skh + inter_cat_neigh + (1|netid), data = data)
  return(model)
}


# cov.delta_maineff = function(n_net, n_memb, rho, p, g, sigma2_y) {
#   #### theoretical variance covariance matrix for delta
#   # rho: ICC of y
#   # g: vector of g_h, length is H
#   # sigma2y: variance of Y
#   c <- 1/(1-rho)
#   d <- -rho/((1-rho)*(1+n_memb*rho))
#   b <- c + d*(n_memb+1)
#   B1 <- rbind(c(p*(c+d), p*b*g),
#               cbind(p*b*g, diag((n_memb+1)*b*g)))
#   B2 <- rbind(n_memb*p*d*g, diag(n_memb*p*b*g))
#   B3 <- t(B2)
#   B4 <- diag(n_memb*p*(c+d*n_memb)*g)
#   return(sigma2_y * solve(B4-B3 %*% solve(B1) %*% B2) /n_net )
# }
# sqrt(diag(cov.delta_maineff(n_net, n_memb, rho=sigma2_u/(sigma2_e+sigma2_u), p, g=g_cat, sigma2_y=sigma2_u+sigma2_e)))
#
# cov.delta = function(n_net, n_memb, rho, p, g, sigma2_y) {
#   #### theoretical variance covariance matrix for delta
#   # rho: ICC of y
#   # g: vector of g_h, length is H
#   # sigma2y: variance of Y
#   c <- 1/(1-rho)
#   d <- -rho/((1-rho)*(1+n_memb*rho))
#   b <- c + d*(n_memb+1)
#   B1 <- diag((n_memb+1)*b*g)
#   B2 <- diag(n_memb*p*b*g)
#   B3 <- t(B2)
#   B4 <- diag(n_memb*p*(c+d*n_memb)*g)
#   return(sigma2_y * solve(B4-B3 %*% solve(B1) %*% B2) /n_net )
# }
# sqrt(diag(cov.delta(n_net, n_memb, rho=sigma2_u/(sigma2_e+sigma2_u), p, g=g_cat, sigma2_y=sigma2_u+sigma2_e)))


### variance too small
## Tukey test
# summary(glht(mod_sim_mixed, mcp(inter_cat_neighb="Tukey")))
# summary(glht(mod_sim_mixed, mcp(inter_cat_neighb="Dunnet")))
# # compare
#
# rep <- 1000
# model_rep <- list()
# true_para <- c(gamma, tau, v_delta)
# bias_est <- sd_para<- matrix(NA, nrow=rep, ncol=length(true_para))
#
# pval_test <- c()
# for (i in 1:rep) {
#   dat <- gen_data(n_net, n_memb, n_cat, g_cat, p, v_zeta, v_delta, sigma2_u, sigma2_e)
#   model <- lmerTest::lmer(y ~ 1 + treat  + net_cat:n_neighb + (1|net_cat), data = dat)
#   est_para <- as.numeric(summary(model)$coef[,"Estimate"])
#   sd_para[i, ] <- as.numeric(summary(model)$coef[,"Std. Error"])
#   bias_est[i, ] <- est_para - true_para
#   pval_test[i] <- anova(model)["net_cat:n_neighb", "Pr(>F)"]
# }
#
# ## report
# names <- c("gamma", "tau", paste0("delta_", 1:4))
# tbl <- data.frame(true_para = true_para, avg_bias = colMeans(bias_est), avg_sd = colMeans(sd_para))
# row.names(tbl) <- names
# # mcb test (confidence intervals)
# # coverage for the best group -- test size (1-alpha)
#
# mod_sim_aov <- aov(y ~ 0 + treat  + net_cat:n_neighb + Error(1/net_cat), data = dat_sim)
# summary(mod_sim_mixed)$coef
#
# tukhsd2mcb(mod_sim_aov,conf.level = 0.95, two.sided = FALSE, method = "MCB")
#

