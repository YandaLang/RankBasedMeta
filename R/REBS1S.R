#' Meta-analysis for random effects - bootstrap
#' @description Meta-analysis for random effects - bootstrap
#' @usage
#' REBS1S(info,alpha=0.05)
#' @param info summary statistics
#' @param iteration number of iterations for bootstrap, default iteration=1000
#' @param alpha significance level, default alpha=0.05
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @note format for summary statistics (point estimate 1, point estimate 2, sample size 1, sample size 2)
#' @examples
#' study1 <- c(247.5, 178.0, 16, 17)
#' study2 <- c(270.0, 167.0, 23, 23)
#' info <- rbind(study1,study2)
#' REBS1S(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export

REBS1S = function(info,iteration=1000,alpha=0.05){
  K = nrow(info)
  est1 = info[,1]; est2 = info[,2]; n1 = info[,3]; n2 = info[,4]
  deltahat = est1-est2

  psq = rep(NA,K)
  for (i in 1:K) {
    psq[i] = (n1[i]*n2[i]/(n1[i]+n2[i]))/sum(n1*n2/(n1+n2))
  }

  est_BS = rep(NA,iteration)
  for (b in 1:iteration){
    ind_BS = sample(1:K,K,replace=TRUE)
    delta_BS = deltahat[ind_BS]
    psq_BS = psq[ind_BS]
    np = sum(psq_BS)
    psq_BS = psq_BS/np
    est_BS[b] = sum(psq_BS*delta_BS)
  }

  deltatilde = mean(est_BS)

  ttvar = var(est_BS)

  lower = qt(alpha/2,K-1)*sqrt(ttvar)+deltatilde
  upper = qt(1-alpha/2,K-1)*sqrt(ttvar)+deltatilde

  combineCI = cbind(deltatilde, lower, upper)
  colnames(combineCI) = c("Estimate","CI.lowerbound","CI.upperbound")
  return(combineCI)
}
