#' Non-Pooled Variance Two-sample Hodges Lehmann estimate
#' @description Calculate two-sample Hodges Lehmann estimate with non-pooled variance
#' @usage
#' get2SHL12(dat, alpha=0.05)
#' @param dat full sample dataframe
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @examples
#' x1 <- rnorm(15)
#' x2 <- rnorm(15)
#' dat <- data.frame(x1,x2)
#' get2SHL12(dat,alpha=0.05)
#' @keywords summary statistics
#' @keywords two-sample Hodges Lehmann estimate

#' @import stats Rfit
#' @export


get2SHL12 = function(dat,alpha=0.05){
  h = ncol(dat);   K = h/2
  delta_hl2samp = rep(NA,K); tausqw = rep(NA,K)
  n_1 = rep(NA,K); n_2 = rep(NA,K)

  # calculate median for ordered diff and its tau
  for(i in 1:K){
    temp1 = dat[,2*i-1]; temp2 = dat[,2*i]
    temp1 = temp1[!is.na(temp1)]; temp2 = temp2[!is.na(temp2)]
    n1 = length(temp1); n2 = length(temp2)

    D = as.vector(outer(temp1,temp2,"-"))
    ds = sort(D)
    delta_hl2samp[i] = median(ds)

    tau1 = gettau(temp1,p=1); tau2 = gettau(temp2,p=1)
    tausqw[i] = tau1^2/n1+tau2^2/n2

    n_1[i] = n1; n_2[i] = n2
  }

  summaryMed = cbind(delta_hl2samp,tausqw,n_1,n_2)
  return(summaryMed)
}
