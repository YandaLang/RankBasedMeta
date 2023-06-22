#' Meta-analysis for fixed effects - LI2
#' @description Meta-analysis for fixed effects - LI2
#' @usage
#' FELI2(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level, default alpha=0.05
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @note format for summary statistics (point estimate 1, point estimate 2)
#' @examples
#' study1 <- c(247.5, 178.0)
#' study2 <- c(270.0, 167.0)
#' study3 <- c(243.0, 186.3)
#' study4 <- c(215.2, 133.7)
#' study5 <- c(245.5, 171.0)
#' info <- rbind(study1,study2,study3,study4,study5)
#' FELI2(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


FELI2 = function(info,alpha=0.05){
  K = nrow(info)
  est1 = info[,1]; est2 = info[,2]
  deltahat = est1-est2

  deltatilde = median(deltahat)

  D = sort(deltahat)
  c = round(K/2-qnorm(1-alpha/2)*sqrt(K/4)-1/2)
  lower = D[c+1]
  upper = D[K-c]

  combineCI = cbind(deltatilde, lower, upper)
  colnames(combineCI) = c("Estimate","CI.lowerbound","CI.upperbound")
  rownames(combineCI) = c(" ")
  return(combineCI)
}
