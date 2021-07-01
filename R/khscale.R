#' Kloke-Yabes procedure
#' @description Kloke-Yabes procedure
#' @usage
#' khscale(estmat,ssmat)
#' @param estmat shift estimates matrix
#' @param ssmat sample size matrix
#' @author Yanda Lang, Joseph McKean
#' @note format for summary statistics (location)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit npsmReg2
#' @export

khscale <- function(estmat,ssmat){
  #
  #   columns of estmate are the estimates xi_{k} of
  #   treatment 1 and 2, rspectively.  ssmat contains
  #   the sample sizes n1_k and n2_k.    So Delta is xi_{1,k} - xi_{2,k}
  #   Scores for now are set at Wilcoxon.
  #

  vht1 <- estmat[,1]; vn1 <- ssmat[,1]
  vht2 <- estmat[,2]; vn2 <- ssmat[,2]

  #  next z is Deltahat
  z <- vht1 - vht2
  x1 <- ((1/vn1) + (1/vn2))^(-.5)
  zstar <- x1*z

  #  This routine is in  library npsmReg2.
  #   Performs a correct fit thru the origin for R regression estimates.
  fit1 <- wtedrb(x1,zstar,scores=Rfit::wscores)
  sighat <- mad(fit1$ehatst)
  sighatk <- sighat/x1
  sighatk2 <- sighatk^2

  xi1 <- sqrt(vn1)/sighat
  xi2 <- sqrt(vn2)/sighat

  # Can't use just shifts input here, but the fix is in the notes K.4
  t1star <- xi1*vht1
  t2star <- xi2*vht2

  yi <- c(); c1i <- c(); c2i <- c(); K <- length(vht1)
  for(i in 1:K){
    c1i <- c(c1i,c(xi1[i],0))
    c2i <- c(c2i,c(0,xi2[i]))
    yi <- c(yi,c(t1star[i],t2star[i]))
  }
  xmati <- cbind(c1i,c2i)
  fit2 <- wtedrb(xmati,yi,scores=Rfit::wscores)
  eistar <- fit2$ehatst

  eistar1 <- rep(0,K); eistar2 <- rep(0,K); are <- rep(0,K)
  for(i in 1:K){
    eistar1[i] <- eistar[2*i-1]/xi1[i]
    eistar2[i] <- eistar[2*i]/xi2[i]
    are[i] <- (eistar1[i] + eistar2[i])/2
  }
  sigahat <- mad(are)
  sigahat2 <- sigahat^2

  wtk <- rep(0,K)
  for(i in 1:K){
    wtk[i] <- 1/(sigahat2 + sighatk2[i])
  }
  #   Output:  sighatk2 is the vector of variances of Delta-hat
  #  i.e  (sighat sqrt((1/n_k^{(1)} + 1/n_k^{(2)}))^2
  #wtk are the wts for standardizing 1/(hat{sigma_a}^2 + sighatk2)
  #  WE DID NOT ESTIMTE tau!!!!!;
  #sigahat2 is hat{sigma_a}^2;   are are the predicted random effects;
  #wtk are the wts for standardizing 1/(hat{sigma_a}^2 + sighatk2)
  #
  list(sighatk2=sighatk2,sigahat2=sigahat2,wtk=wtk,are=are)
}
