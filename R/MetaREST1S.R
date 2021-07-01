#' Meta-analysis for random-effect - traditional
#' @description Meta-analysis for random-effect - traditional
#' @usage
#' MetaREST1S(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @note format for summary statistics (location, scale, sample size)
#' @examples
#' study1_treat <- c(178.0, 109.5, 17)
#' study1_cont <- c(247.5, 113.8, 16)
#' study2_treat <- c(167.0, 162.6, 23)
#' study2_cont <- c(270.0, 93.8, 23)
#' info <- rbind(study1_treat,study1_cont,study2_treat,study2_cont)
#' MetaREST1S(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


MetaREST1S = function(info,alpha=0.05){
  # Define h studies split into K pair of studies
  h = nrow(info)
  K = h/2
  # Define est, tauhat, n, alpha
  est = info[,1]; tauhat = info[,2]; n = info[,3]
  # Define estimator delta hat - deltahat, and its variance - tauhatsq
  deltahat = rep(NA,K)
  tauhatsq = rep(NA,K)

  # Calculate shift between two sample median, and its variance
  for (i in 1:K) {
    j = 2*i-1
    # shift
    deltahat[i] = est[j]-est[j+1]
    # variance
    tauhatsq[i] = tauhat[j]^2/n[j]+tauhat[j+1]^2/n[j+1]
  }

  weightk = 1/tauhatsq
  # The estimate of delta based on the weighted average of the K summary, deltatilde
  deltatilde = sum(deltahat/tauhatsq)/sum(1/tauhatsq)

  # Total variance (between+within)
  Qsq = sum((deltahat-deltatilde)^2*weightk)

  # Define coefficient c
  c = sum(weightk)-sum(weightk^2)/sum(weightk)

  # Calculate variance for between
  if (Qsq > K-1){
    tauhatsq_between = (Qsq-(K-1))/c
  } else {
    tauhatsq_between = 0
  }

  # New weight
  wstar = 1/(tauhatsq+tauhatsq_between)

  # Random effect estimate
  deltatildeR = sum(wstar*deltahat)/sum(wstar)

  # Jackknife
  # Define estimate of delta based on K-1 combined summary statistics, deltatilde_kminus1;
  # Define Jackknife replicates for deltatilde, deltatilde_knife
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

  # Jackknife variance estimate of deltatildeR, B_hat
  # mean of Jackknife replicates for deltatildeR, mean_deltatildeR_knife
  mean_deltatildeR_knife = sum(deltatildeR_knife)/K
  B_hat = sum((deltatildeR_knife-mean_deltatildeR_knife)^2)/(K*(K-1))

  # proposed combined CI
  lower = qt(alpha/2,K-1)*sqrt(B_hat)+deltatildeR
  upper = qt(1-alpha/2,K-1)*sqrt(B_hat)+deltatildeR

  # collect values
  sqrtB = sqrt(B_hat)
  combineCI = c(lower, upper)
  return(combineCI)
}
