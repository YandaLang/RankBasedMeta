#' Meta-analysis for random-effect - D4
#' @description Meta-analysis for random-effect - D4
#' @usage
#' MetaRED4(info,alpha=0.05)
#' @param info Hodges Lehmann estimates, standard deviation, and sampel size
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @note format for info (HL,sd,sample size)
#' @examples
#' study1_treat <- c(178.0,198.9,17)
#' study1_cont <- c(247.5,179.2,16)
#' study2_treat <- c(167.0,118.8,23)
#' study2_cont <- c(270.0,90.8,23)
#' info <- rbind(study1_treat,study1_cont,study2_treat,study2_cont)
#' MetaRED4(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


MetaRED4 = function(info,alpha=0.05){
  # Define h studies split into K pair of studies
  h = nrow(info)
  K = h/2
  # Define median, n, alpha
  est = info[,1]; sigmahat = info[,2]; n = info[,3]
  # Define estimator delta hat - deltahat, and its variance - tauhatsq
  deltahat = rep(NA,K)
  sigmahatk = rep(NA,K)

  # Calculate shift between two sample median, and its variance
  for (i in 1:K) {
    j = 2*i-1
    # shift
    deltahat[i] = est[j]-est[j+1]
    # variance
    sigmahatk[i] = sqrt(sigmahat[j]^2/n[j]+sigmahat[j+1]^2/n[j+1])
  }
  tauhatk = sigmahatk * sqrt(3/pi)

  weightk = 1/tauhatk^2
  # The estimate of delta based on the weighted average of the K summary, deltatilde
  deltatilde = sum(deltahat/tauhatk^2)/sum(1/tauhatk^2)

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

  tauhatsq_within = tauhatk^2
  # New weight
  wstar = 1/(tauhatsq_within+tauhatsq_between)

  # Random effect estimate
  deltatildeR = sum(wstar*deltahat)/sum(wstar)
  # intraclass correction coefficient estimate
  # ICCE = tauhatsq_between/(tauhatsq_within+tauhatsq_between)

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

  sqrtB = sqrt(B_hat)
  combineCI = c(lower, upper)
  return(combineCI)
}
