#' Pooled Variance Two-sample Hodges Lehmann estimate with HLPT scale parameter
#' @description Calculate two-sample Hodges Lehmann estimate with pooled variance with HLPT scale parameter
#' @usage
#' get2SHL12P(dat, alpha=0.05)
#' @param dat full sample dataframe
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @examples
#' x1 <- rnorm(15)
#' x2 <- rnorm(15)
#' dat <- data.frame(x1,x2)
#' get2SHL12P(dat,alpha=0.05)
#' @keywords summary statistics
#' @keywords two-sample Hodges Lehmann estimate

#' @import stats Rfit
#' @export


get2SHL12P = function(dat,alpha=0.05){
  h = ncol(dat);   K = h/2
  delta_hl2samp = rep(NA,K); tausqw = rep(NA,K)
  n_1 = rep(NA,K); n_2 = rep(NA,K)

  # calculate median for ordered diff and its tau
  for(i in 1:K){
    temp1 = dat[,2*i-1]; temp2 = dat[,2*i]
    temp1 = temp1[!is.na(temp1)]; temp2 = temp2[!is.na(temp2)]
    n1 = length(temp1); n2 = length(temp2); n = n1+n2

    D = as.vector(outer(temp1,temp2,"-"))
    ds = sort(D)
    delta_hl2samp[i] = median(ds)

    temp_walsh1 = sort(walsh(temp1))
    k_walsh1 = round(n1*(n1+1)/4-0.5-qnorm(0.975)*sqrt(n1*(n1+1)*(2*n1+1)/24))
    CI_lower1 = temp_walsh1[k_walsh1+1]
    CI_upper1 = temp_walsh1[n1*(n1+1)/2-k_walsh1]
    tauhat1 = sqrt(n1)*(CI_upper1-CI_lower1)/(2*abs(qnorm(alpha/2)))

    temp_walsh2 = sort(walsh(temp2))
    k_walsh2 = round(n2*(n2+1)/4-0.5-qnorm(0.975)*sqrt(n2*(n2+1)*(2*n2+1)/24))
    CI_lower2 = temp_walsh2[k_walsh2+1]
    CI_upper2 = temp_walsh2[n2*(n2+1)/2-k_walsh2]
    tauhat2 = sqrt(n2)*(CI_upper2-CI_lower2)/(2*abs(qnorm(alpha/2)))

    tauhatsq = ((n1-1)*tauhat1^2+(n2-1)*tauhat2^2)/(n1+n2-2)
    tausqw[i] = tauhatsq*(1/n1+1/n2)

    n_1[i] = n1; n_2[i] = n2
  }

  summaryMed = cbind(delta_hl2samp,tausqw,n_1,n_2)
  return(summaryMed)
}

