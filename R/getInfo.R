#' Summary statistics
#' @description Calculate summary statistics
#' @usage
#' getInfo(dat, alpha=0.05)
#' @param dat full sample dataframe
#' @param alpha significance level for tau estimation
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
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
  h = ncol(dat); K=h/2
  xbar = rep(NA,h); stdev = rep(NA,h)
  med = rep(NA,h); tauhat_mm = rep(NA,h)
  med_walsh = rep(NA,h); tauhat_KSM = rep(NA,h)
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
    med_walsh[i] = median(temp_walsh)
    tauhat_KSM[i] = gettau(temp,p=1)
  }

  n1 = rep(NA,K); n2 = rep(NA,K)
  xbar1 = rep(NA,K); xbar2 = rep(NA,K)
  stdev1 = rep(NA,K); stdev2 = rep(NA,K)
  med1 = rep(NA,K); med2 = rep(NA,K)
  tauhat_mm1 = rep(NA,K); tauhat_mm2 = rep(NA,K)
  med_walsh1 = rep(NA,K); med_walsh2 = rep(NA,K)
  tauhat_KSM1 = rep(NA,K); tauhat_KSM2 = rep(NA,K)
  for(i in 1:K){
    n1[i] = n[2*i-1]; n2[i] = n[2*i]
    xbar1[i] = xbar[2*i-1]; xbar2[i] = xbar[2*i]
    stdev1[i] = stdev[2*i-1]; stdev2[i] = stdev[2*i]
    med1[i] = med[2*i-1]; med2[i] = med[2*i]
    tauhat_mm1[i] = tauhat_mm[2*i-1]; tauhat_mm2[i] = tauhat_mm[2*i]
    med_walsh1[i] = med_walsh[2*i-1]; med_walsh2[i] = med_walsh[2*i]
    tauhat_KSM1[i] = tauhat_KSM[2*i-1]; tauhat_KSM2[i] = tauhat_KSM[2*i]
  }

  delta_w = rep(NA,K); tauhat_wk = rep(NA,K); tauhat_wp = rep(NA,K)
  n=n1+n2

  for(i in 1:K){
    temp1 = dat[,2*i-1]; temp2 = dat[,2*i]
    temp1 = temp1[!is.na(temp1)]; temp2 = temp2[!is.na(temp2)]

    tempz = c(temp1,temp2)
    ind = c(rep(0,n1[i]),rep(1,n2[i]))
    fit = rfit(tempz~ind)
    delta_w[i] = fit$coef[2]
    tauhat_wk[i] = fit$tauhat

    tauhat_wp[i] = sqrt(tauhat_wk[i]^2*(1/n1[i]+1/n2[i]))
  }

  summaryMed = cbind(n1,n2,xbar1,xbar2,stdev1,stdev2,med1,med2,tauhat_mm1,tauhat_mm2,med_walsh1,med_walsh2,tauhat_KSM1,tauhat_KSM2,delta_w,tauhat_wp)
  return(summaryMed)
}



