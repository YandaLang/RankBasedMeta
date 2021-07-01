#' Meta-analysis for fixed-effect non-pooled variance
#' @description Meta-analysis for fixed-effect non-pooled variance
#' @usage
#' MetaFE1S2T(info,alpha=0.05)
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
#' MetaFE1S2T(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


MetaFE1S2T = function(info,alpha=0.05){
  # Define h studies split into K pair of studies
  h = nrow(info)
  K = h/2
  # Define median, tauhat, n
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

