#' Meta-analysis for fixed effects - two-sample Wilcoxon non-pooled variance
#' @description Meta-analysis for fixed effects - two-sample Wilcoxon non-pooled variance
#' @usage
#' FEWNP(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @note format for summary statistics (shift, scale 1, scale 2, sample size 1, sample size 2)
#' @examples
#' study1 <- c(69.5, 19.5, 13.8, 17, 16)
#' study2 <- c(103, 12.6, 9.8, 20, 30)
#' info <- rbind(study1,study2)
#' FEWNP(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export

FEWNP = function(info,alpha=0.05){
  K = nrow(info)
  deltahat = info[,1]; tauhat1 = info[,2]; tauhat2 = info[,3]; n1 = info[,4]; n2 = info[,5]
  tauhatsq1 = tauhat1^2; tauhatsq2 = tauhat2^2
  tauhatsq = tauhatsq1/n1+tauhatsq2/n2

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
  colnames(combineCI) = c("DeltaEstimate","LowerBound","UpperBound")
  return(combineCI)
}

