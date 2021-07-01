#' Meta-analysis for random-effect - two-sample KY
#' @description Meta-analysis for random-effect - two-sample KY
#' @usage
#' MetaREKY2S(info,alpha=0.05)
#' @param info shift estimates, scale estimates, location parameters, and sampel size
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @note format for info (location,scale,location parameter1,location parameter2,sample size 1,sample size 2)
#' @examples
#' study1_2SHL <- c(-80.0, 1243.9, 178.0, 247.5, 17, 16)
#' study2_2SHL <- c(-89.0, 1029.3, 167.0, 270.0, 23, 23)
#' info <- rbind(study1_2SHL,study2_2SHL)
#' MetaREKY2S(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit npsmReg2
#' @export

MetaREKY2S = function(info,alpha=0.05){
  # Define h studies split into K pair of studies
  K = nrow(info)
  # Define shift estimate, sample size
  deltahat = info[,1]; ssmat = info[,5:6]; mumat = info[,3:4]; tauhatsq = info[,2]

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

  # Kloke estimate random effect variances
  fit = khscale2S(deltahat,mumat,ssmat)

  # collect hat{sigma}_a^2
  tauhatsq_within = fit$sighatk2

  # New weight
  wstar = 1/(tauhatsq_within+tauhatsq_between)

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
