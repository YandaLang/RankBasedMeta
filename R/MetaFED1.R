#' Meta-analysis for fixed-effect - D1
#' @description Meta-analysis for fixed-effect - D1
#' @usage
#' MetaFED1(info,alpha=0.05)
#' @param info Hodges Lehmann estimates
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @note format for info (Hodges Lehmann estimates)
#' @examples
#' study1_treat <- c(178.0)
#' study1_cont <- c(247.5)
#' study2_treat <- c(167.0)
#' study2_cont <- c(270.0)
#' info <- rbind(study1_treat,study1_cont,study2_treat,study2_cont)
#' MetaFED1(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export

MetaFED1 = function(info,alpha=0.05){
  # Define h studies split into k pair of studies
  K = nrow(info)/2
  # Define median
  med = info[,1]
  # Define estimator delta hat - deltahat
  deltahat = rep(NA,K)

  # Calculate shift between two sample median
  for (i in 1:K) {
    j = 2*i-1
    # shift
    deltahat[i] = med[j]-med[j+1]
  }

  # The estimate of delta using median of differences
  deltatilde = median(deltahat)

  H = sort(deltahat)
  c = round(K/2-qnorm(1-alpha/2)*sqrt(K/4)-1/2)

  # proposed combined CI
  lower = H[c+1]
  upper = H[K-c]

  # collect values
  combineCI = c(lower, upper)
  return(combineCI)
}
