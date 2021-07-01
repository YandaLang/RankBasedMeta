#' Meta-analysis for fixed-effect - C1
#' @description Meta-analysis for fixed-effect - C1
#' @usage
#' MetaFEC1(info,alpha=0.05)
#' @param info sample median
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @note format for info (median)
#' @examples
#' study1_treat <- c(178.0)
#' study1_cont <- c(247.5)
#' study2_treat <- c(167.0)
#' study2_cont <- c(270.0)
#' info <- rbind(study1_treat,study1_cont,study2_treat,study2_cont)
#' MetaFEC1(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


MetaFEC1 = function(info,alpha=0.05){
  # Define h studies split into k pair of studies
  h = nrow(info)
  K = h/2
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

  D = sort(deltahat)
  c = round(K/2-qnorm(1-alpha/2)*sqrt(K/4)-1/2)

  # proposed combined CI
  lower = D[c+1]
  upper = D[K-c]

  # collect values
  combineCI = c(lower, upper)
  return(combineCI)
}
