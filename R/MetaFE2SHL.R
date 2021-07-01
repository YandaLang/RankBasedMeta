#' Meta-analysis for fixed-effect - two-sample Hodges Lehmann
#' @description Meta-analysis for fixed-effect - two-sample Hodges Lehmann
#' @usage
#' MetaFE2SHL(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @note format for summary statistics (location, scale)
#' @examples
#' study1 <- c(-80.0, 1243.9)
#' study2 <- c(-89.0, 1029.3)
#' info <- rbind(study1,study2)
#' MetaFE2SHL(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export

MetaFE2SHL = function(info,alpha=0.05){
  # Define h studies split into K pair of studies
  K = nrow(info)
  # Define median diff, tauhatsq
  deltahat = info[,1]; tauhatsq = info[,2]

  # Jackknife
  # The estimate of delta based on the weighted average of the K summary, deltatilde
  deltatilde = sum(deltahat/tauhatsq)/sum(1/tauhatsq)

  # Define estimate of delta based on K-1 combined summary statistics, deltatilde_kminus1;
  # Define Jackknife replicates for deltatilde, deltatilde_knife
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

  # Jackknife variance estimate of deltatilde, B_hat
  # mean of Jackknife replicates for deltatilde, mean_deltatilde_knife
  mean_deltatilde_knife = sum(deltatilde_knife)/K
  B_hat = sum((deltatilde_knife-mean_deltatilde_knife)^2)/(K*(K-1))

  # proposed combined CI
  lower = qt(alpha/2,K-1)*sqrt(B_hat)+deltatilde
  upper = qt(1-alpha/2,K-1)*sqrt(B_hat)+deltatilde

  sqrtB = sqrt(B_hat)
  combineCI = c(lower, upper)
  return(combineCI)
}
