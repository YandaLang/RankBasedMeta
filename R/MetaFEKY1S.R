#' Meta-analysis for fixed-effect - KY
#' @description Meta-analysis for fixed-effect - KY
#' @usage
#' MetaFEKY1S(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @note format for summary statistics (location,sample size)
#' @examples
#' study1_treat <- c(178.0,17)
#' study1_cont <- c(247.5,16)
#' study2_treat <- c(167.0,23)
#' study2_cont <- c(270.0,23)
#' info <- rbind(study1_treat,study1_cont,study2_treat,study2_cont)
#' MetaFEKY1S(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit npsmReg2
#' @export

MetaFEKY1S = function(info,alpha=0.05){
  # Define h studies split into K pair of studies
  h = nrow(info)
  K = h/2
  # Define median, n
  med = info[,1]; n = info[,2]

  # define estimate, and its sample size
  estmat = matrix(NA,nrow = K,ncol=2)
  nmat = matrix(NA,nrow = K,ncol=2)
  deltahat = rep(NA,K)

  for (i in 1:K) {
    j = 2*i-1
    # est
    estmat[i,1] = med[j]
    estmat[i,2] = med[j+1]
    # sample size
    nmat[i,1] = n[j]
    nmat[i,2] = n[j+1]
    # shift
    deltahat[i] = med[j]-med[j+1]
  }

  # Kloke estimate random effect variances
  fit = khscale(estmat,nmat)

  # New weight
  wstar = fit$sighatk2

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
