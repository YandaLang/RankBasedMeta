#' Meta-analysis for random effects - two-sample Wilcoxon non-pooled variance
#' @description Meta-analysis for random effects - two-sample Wilcoxon non-pooled variance
#' @usage
#' RESTWNP(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @note format for summary statistics (shift, scale 1, scale 2, sample size 1, sample size 2)
#' @examples
#' study1 <- c(69.5, 19.5, 13.8, 17, 16)
#' study2 <- c(103, 12.6, 9.8, 20, 30)
#' info <- rbind(study1,study2)
#' RESTWNP(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


RESTWNP = function(info,alpha=0.05){
  K = nrow(info)
  deltahat = info[,1]; tauhat1 = info[,2]; tauhat2 = info[,3]; n1 = info[,4]; n2 = info[,5]
  tauhatsq1 = tauhat1^2; tauhatsq2 = tauhat2^2
  tauhatsq = tauhatsq1/n1+tauhatsq2/n2

  deltatilde = sum(deltahat/tauhatsq)/sum(1/tauhatsq)

  Qsq = sum((deltahat-deltatilde)^2/tauhatsq)

  c = sum(1/tauhatsq)-sum((1/tauhatsq)^2)/sum(1/tauhatsq)

  if (Qsq > K-1){
    tauhatsq_between = (Qsq-(K-1))/c
  } else {
    tauhatsq_between = 0
  }

  wstar = 1/(tauhatsq+tauhatsq_between)

  deltatildeR = sum(wstar*deltahat)/sum(wstar)

  deltatildeR_kminus1 = rep(NA,K)
  deltatildeR_knife = rep(NA,K)

  for (i in 1:K) {
    deltaknife = deltahat[-i]
    wstarknife = wstar[-i]
    deltatildeR_kminus1[i] = sum(deltaknife*wstarknife)/sum(wstarknife)
  }

  for (i in 1:K) {
    deltatildeR_knife[i] = K*deltatildeR-(K-1)*deltatildeR_kminus1[i]
  }

  mean_deltatildeR_knife = sum(deltatildeR_knife)/K
  B_hat = sum((deltatildeR_knife-mean_deltatildeR_knife)^2)/(K*(K-1))

  lower = qt(alpha/2,K-1)*sqrt(B_hat)+deltatildeR
  upper = qt(1-alpha/2,K-1)*sqrt(B_hat)+deltatildeR

  sqrtB = sqrt(B_hat)
  combineCI = cbind(deltatilde, lower, upper)
  colnames(combineCI) = c("DeltaEstimate","LowerBound","UpperBound")
  return(combineCI)
}

