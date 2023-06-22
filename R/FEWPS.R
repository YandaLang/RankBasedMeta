#' Meta-analysis for fixed effects - two-sample Wilcoxon pooled
#' @description Meta-analysis for fixed effects - two-sample Wilcoxon pooled
#' @usage
#' FEWPS(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level, default alpha=0.05
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @note format for summary statistics (shift, scale)
#' @examples
#' study1 <- c(69.5, 11.65)
#' study2 <- c(103, 18.2)
#' info <- rbind(study1,study2)
#' FEWPS(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export

FEWPS = function(info,alpha=0.05){
  K = nrow(info)
  deltahat = info[,1]; tauhat = info[,2]
  tauhatsq = tauhat^2

  deltatilde = sum(deltahat/tauhatsq)/sum(1/tauhatsq)

  deltatilde_kminus1 = rep(NA,K)
  deltatilde_knife = rep(NA,K)

  for (i in 1:K) {
    deltaknife = deltahat[-i]
    tausqknife = tauhatsq[-i]
    deltatilde_kminus1[i] = sum(deltaknife/tausqknife)/sum(1/tausqknife)
  }

  for (i in 1:K) {
    deltatilde_knife[i] = K*deltatilde-(K-1)*deltatilde_kminus1[i]
  }
  mean_deltatilde_knife = sum(deltatilde_knife)/K
  B_hat = sum((deltatilde_knife-mean_deltatilde_knife)^2)/(K*(K-1))

  lower = qt(alpha/2,K-1)*sqrt(B_hat)+deltatilde
  upper = qt(1-alpha/2,K-1)*sqrt(B_hat)+deltatilde

  sqrtB = sqrt(B_hat)
  combineCI = cbind(deltatilde, lower, upper)
  colnames(combineCI) = c("Estimate","CI.lowerbound","CI.upperbound")
  return(combineCI)
}

