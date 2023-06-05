#' Meta-analysis for fixed effects - LI
#' @description Meta-analysis for fixed effects - LI
#' @usage
#' FELI(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level, default alpha=0.05
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @note format for summary statistics (point estimate 1, point estimate 2, sample size 1, sample size 2)
#' @examples
#' study1 <- c(247.5, 178.0, 16, 17)
#' study2 <- c(270.0, 167.0, 23, 23)
#' info <- rbind(study1,study2)
#' FELI(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


FELI = function(info,alpha=0.05){
  K = nrow(info)
  est1 = info[,1]; est2 = info[,2]; n1 = info[,3]; n2 = info[,4]
  deltahat = est1-est2

  m = rep(NA,K)
  for (i in 1:K) {
    m[i] = (n1[i]*n2[i]/(n1[i]+n2[i]))/sum(n1*n2/(n1+n2))
  }
  deltatilde = sum(m*deltahat)

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

  mean_deltatilde_knife = sum(deltatilde_knife)/K
  B_hat = sum((deltatilde_knife-mean_deltatilde_knife)^2)/(K*(K-1))

  lower = qnorm(alpha/2)*sqrt(B_hat)+deltatilde
  upper = qnorm(1-alpha/2)*sqrt(B_hat)+deltatilde

  sqrtB = sqrt(B_hat)
  combineCI = cbind(deltatilde, lower, upper)
  colnames(combineCI) = c("DeltaEstimate","LowerBound","UpperBound")
  return(combineCI)
}
