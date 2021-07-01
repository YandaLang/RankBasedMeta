#' Summary statistics
#' @description Calculate summary statistics
#' @usage
#' getInfo(dat, alpha=0.05)
#' @param dat full sample dataframe
#' @param alpha significance level
#' @author Yanda Lang, Joseph McKean
#' @examples
#' x1 <- rnorm(15)
#' x2 <- rnorm(15)
#' dat <- data.frame(x1,x2)
#' getInfo(dat,alpha=0.05)
#' @keywords summary statistics
#' @keywords Hodges Lehmann estimate

#' @import stats Rfit
#' @export

getInfo = function(dat,alpha=0.05){
  h = ncol(dat)
  xbar = rep(NA,h); stdev = rep(NA,h)
  med = rep(NA,h); tauhat_mm = rep(NA,h)
  med_walsh1 = rep(NA,h); tauhat_walsh1 = rep(NA,h)
  k_walsh = rep(NA,h); CI_lower = rep(NA,h); CI_upper = rep(NA,h)
  med_walsh2 = rep(NA,h); tauhat_walsh2 = rep(NA,h)
  n = rep(NA,h)

  for(i in 1:h){
    temp = dat[,i]
    temp = temp[!is.na(temp)]
    n[i] = length(temp)
    # LS
    xbar[i] = mean(temp)
    stdev[i] = sd(temp)
    # MM
    med[i] = median(temp)
    tauhat_mm[i] = taustar(temp,0)
    # HLKSM
    temp_walsh = sort(walsh(temp))
    med_walsh1[i] = median(temp_walsh)
    tauhat_walsh1[i] = gettau(temp,p=1)
    # HLCI
    k_walsh[i] = round(n[i]*(n[i]+1)/4-0.5-qnorm(0.975)*sqrt(n[i]*(n[i]+1)*(2*n[i]+1)/24))
    CI_lower[i] = temp_walsh[k_walsh[i]+1]
    CI_upper[i] = temp_walsh[n[i]*(n[i]+1)/2-k_walsh[i]]
    med_walsh2[i] = (CI_lower[i]+CI_upper[i])/2
    tauhat_walsh2[i] = sqrt(n[i])*(CI_upper[i]-CI_lower[i])/(2*abs(qnorm(alpha/2)))
  }
  summaryMed = cbind(alpha,n,xbar,stdev,med,tauhat_mm,med_walsh1,tauhat_walsh1,
                     med_walsh2,tauhat_walsh2)
  return(summaryMed)
}


