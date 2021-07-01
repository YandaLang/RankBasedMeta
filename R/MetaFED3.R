#' Meta-analysis for fixed-effect - D3
#' @description Meta-analysis for fixed-effect - D3
#' @usage
#' MetaFED3(info,alpha=0.05)
#' @param info Hodges Lehmann estimates
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @note format for info (Hodges Lehmann estimates, sample size)
#' @examples
#' study1_treat <- c(178.0, 17)
#' study1_cont <- c(247.5, 16)
#' study2_treat <- c(167.0, 23)
#' study2_cont <- c(270.0, 23)
#' info <- rbind(study1_treat,study1_cont,study2_treat,study2_cont)
#' MetaFED3(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


MetaFED3 = function(info,alpha=0.05){
  # Define h studies split into k pair of studies
  h = nrow(info)
  K = h/2
  # Define median, n
  med = info[,1]; n = info[,2]
  # Define estimator delta hat - deltahat, sample size, n1, n2
  deltahat = rep(NA,K)
  n1 = rep(NA,K); n2 = rep(NA,K)

  # Calculate shift between two sample median, and sample size
  for (i in 1:K) {
    j = 2*i-1
    # shift
    deltahat[i] = med[j]-med[j+1]
    # sample size n1, n2
    n1[i] = n[j]; n2[i] = n[j+1]
  }

  # Jackknife
  # The estimate of delta based on the median of the weighted differences
  m = rep(NA,K)
  for (i in 1:K) {
    m[i] = sqrt(n1[i]*n2[i]/(n1[i]+n2[i]))/sum(sqrt(n1*n2/(n1+n2)))
  }
  deltatilde = sum(m*deltahat)

  # Define estimate of delta based on K-1 combined summary statistics, deltatilde_kminus1;
  # Define Jackknife replicates for deltatilde, deltatilde_knife
  deltatilde_kminus1 = rep(NA,K)
  deltatilde_knife = rep(NA,K)
  mknife = rep(NA,K-1)

  for (i in 1:K) {
    deltaknife = deltahat[-i]
    n1knife = n1[-i]
    n2knife = n2[-i]
    for (j in 1:K-1) {
      mknife[j] = sqrt(n1knife[j]*n2knife[j]/(n1knife[j]+n2knife[j]))/sum(sqrt(n1knife*n2knife/(n1knife+n2knife)))
    }
    deltatilde_kminus1[i] = sum(mknife*deltaknife)
  }

  for (i in 1:K) {
    deltatilde_knife[i] = K*deltatilde-(K-1)*deltatilde_kminus1[i]
  }

  # Jackknife variance estimate of deltatilde, B_hat
  # mean of Jackknife replicates for deltatilde, mean_deltatilde_knife
  mean_deltatilde_knife = sum(deltatilde_knife)/K
  B_hat = sum((deltatilde_knife-mean_deltatilde_knife)^2)/(K*(K-1))

  # proposed combined CI
  lower = qnorm(alpha/2)*sqrt(B_hat)+deltatilde
  upper = qnorm(1-alpha/2)*sqrt(B_hat)+deltatilde

  # collect values
  sqrtB = sqrt(B_hat)
  combineCI = c(lower, upper)
  return(combineCI)
}
