#' Meta-analysis for fixed effects pooled variance
#' @description Meta-analysis for fixed effects pooled variance
#' @usage
#' FE1SPT(info,alpha=0.05)
#' @param info summary statistics
#' @param alpha significance level, default alpha=0.05
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @note format for summary statistics (point estimate 1, point estimate 2, scale 1, scale 2, sample size 1, sample size 2)
#' @examples
#' study1 <- c(247.5, 178.0, 13.8, 19.5, 16, 17)
#' study2 <- c(270.0, 167.0, 9.8, 12.6, 23, 23)
#' info <- rbind(study1,study2)
#' FE1SPT(info,alpha=0.05)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export


FE1SPT = function(info,alpha=0.05){
  K = nrow(info)
  est1 = info[,1]; est2 = info[,2]; tauhat1 = info[,3]; tauhat2 = info[,4]; n1 = info[,5]; n2 = info[,6]

  deltahat = est1-est2
  taup = ((n1-1)*tauhat1^2+(n2-1)*tauhat2^2)/(n1+n2-2)
  tauhatsq = taup*sqrt(1/n1+1/n2)

  deltatilde = sum(deltahat/tauhatsq)/sum(1/tauhatsq)

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

  mean_deltatilde_knife = sum(deltatilde_knife)/K
  B_hat = sum((deltatilde_knife-mean_deltatilde_knife)^2)/(K*(K-1))

  lower = qt(alpha/2,K-1)*sqrt(B_hat)+deltatilde
  upper = qt(1-alpha/2,K-1)*sqrt(B_hat)+deltatilde

  sqrtB = sqrt(B_hat)
  combineCI = cbind(deltatilde, lower, upper)
  colnames(combineCI) = c("Estimate","CI.lowerbound","CI.upperbound")
  return(combineCI)
}

