#' Meta-analysis for fixed-effect - E1
#' @description Meta-analysis for fixed-effect - E1
#' @usage
#' MetaFEE1(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @note format for info (2SHL,sample size 1,sample size 2)
#' @examples
#' study1 <- c(-80.0, 17, 16)
#' study2 <- c(-89.0, 23, 23)
#' info <- rbind(study1,study2)
#' MetaFEE1(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


MetaFEE1 = function(info,alpha=0.05){
  # Define h studies split into K pair of studies
  K = nrow(info)
  # Define median diff, sample size
  deltahat = info[,1]; n1 = info[,2]; n2 = info[,3]

  # The estimate of delta
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
