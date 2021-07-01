#' Meta-analysis for random-effect - two-sample traditional
#' @description Meta-analysis for random-effect - two-sample traditional
#' @usage
#' MetaREST2S(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @note format for summary statistics (location, scale, sample size)
#' @examples
#' study1 <- c(-80.0, 1243.9)
#' study2 <- c(-89.0, 1029.3)
#' info <- rbind(study1,study2)
#' MetaREST2S(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export

MetaREST2S = function(info,alpha=0.05){
  # Define h studies split into k pair of studies
  K = nrow(info)
  # Define median diff, tauhatsq
  deltahat = info[,1]; tauhatsq = info[,2]

  # The estimate of delta based on the weighted average of the K summary, deltatilde
  deltatilde = sum(deltahat/tauhatsq)/sum(1/tauhatsq)

  # Total variance (between+within)
  Qsq = sum((deltahat-deltatilde)^2/tauhatsq)

  # Define coefficient c
  c = sum(1/tauhatsq)-sum((1/tauhatsq)^2)/sum(1/tauhatsq)

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

