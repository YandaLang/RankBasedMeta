#' Meta-analysis q-q Plots
#' @description q-q plots based on selected residuals
#' @usage
#' meta.qqplot(data, comparison = FALSE, effect = "fixed", residual = "MMPS")
#' @param data input data
#' @param comparison indicating q-q plots for comparison over all studies or one q-q plot for the entire study; default comparison q-q plots
#' @param effect fixed effects or random effects
#' @param residual Meta-analysis methods, default MMNP
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @examples
#' sample1_treat <- rnorm(25,0,1)
#' sample1_con <- rnorm(25,0,1)
#' sample2_treat <- rnorm(25,1,1)
#' sample2_con <- rnorm(25,1,1)
#' dat <- data.frame(sample1_treat,sample1_con,sample2_treat,sample2_con)
#' meta.qqplot(dat, effect = "fixed", residual = "MMPS")
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export

meta.qqplot = function(data, comparison = TRUE, effect = "fixed", residual = "MMPS"){
  dat_col = ncol(data); K = dat_col/2
  Yhat = rmeta(data, effect = effect, method = residual)[[1]][1,1]
  if (comparison == TRUE){
    par(mfrow=c(K,2))
    for (i in 1:dat_col){
      Yresidual = data[,i]-Yhat
      qqnorm(Yresidual,main=residual);qqline(Yresidual)
    }
    par(mfrow=c(1,1))
  }else if (comparison == FALSE){
    Y = as.vector(t(data))
    Y=Y[!is.na(Y)]
    Yresidual = Y-Yhat
    qqnorm(Yresidual,main=residual);qqline(Yresidual)
  }
}



