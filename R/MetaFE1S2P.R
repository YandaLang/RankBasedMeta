#' Meta-analysis for fixed-effect pooled variance
#' @description Meta-analysis for fixed-effect pooled variance
#' @usage
#' MetaFE1SPT(info,alpha=0.05)
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
#' MetaFE1SPT(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


MetaFE1SPT = function(info,alpha=0.05){
  # Define h studies split into k pair of studies
  h = nrow(info)
  k = h/2
  # Define median, tauhat, n, alpha
  est = info[,1]; tauhat = info[,2]; n = info[,3]
  # Define estimator delta hat - deltahat, and its pooled variance - tauhatsq
  deltahat = rep(NA,k)
  taup = rep(NA,k)
  tauhatsq = rep(NA,k)

  # Calculate shift between two sample median, and its pooled variance
  for (i in 1:k) {
    j = 2*i-1
    # shift
    deltahat[i] = est[j]-est[j+1]
    # pooled variance
    taup[i] = ((n[j]-1)*tauhat[j]^2+(n[j+1]-1)*tauhat[j+1]^2)/(n[j]+n[j+1]-2)
    tauhatsq[i] = taup[i]*sqrt(1/n[j]+1/n[j+1])
  }

  # Jackknife
  # The estimate of delta based on the weighted average of the K summary, deltatilde
  deltatilde = sum(deltahat/tauhatsq)/sum(1/tauhatsq)

  # Define estimate of delta based on K-1 combined summary statistics, deltatilde_kminus1;
  # Define Jackknife replicates for deltatilde, deltatilde_knife
  deltatilde_kminus1 = rep(NA,k)
  deltatilde_knife = rep(NA,k)

  for (i in 1:k) {
    deltaknife = deltahat[-i]
    tausqknife = tauhatsq[-i]
    deltatilde_kminus1[i] = sum(deltaknife/tausqknife)/sum(1/tausqknife)
  }

  for (i in 1:k) {
    deltatilde_knife[i] = k*deltatilde-(k-1)*deltatilde_kminus1[i]
  }

  # Jackknife variance estimate of deltatilde, B_hat
  # mean of Jackknife replicates for deltatilde, mean_deltatilde_knife
  mean_deltatilde_knife = sum(deltatilde_knife)/k
  B_hat = sum((deltatilde_knife-mean_deltatilde_knife)^2)/(k*(k-1))

  # proposed combined CI
  lower = qt(alpha/2,k-1)*sqrt(B_hat)+deltatilde
  upper = qt(1-alpha/2,k-1)*sqrt(B_hat)+deltatilde

  sqrtB = sqrt(B_hat)
  combineCI = c(lower, upper)
  return(combineCI)
}

