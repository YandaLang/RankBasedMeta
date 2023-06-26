#' Rank-based Meta-analysis
#' @description Rank-based Meta-analysis
#' @usage
#' rmeta(data, effect = "fixed", method = "MMNP", fullsample = TRUE)
#' @param data input data
#' @param effect fixed effects or random effects
#' @param method Meta-analysis methods, default MMNP
#' @param fullsample data type. full sample or summary statistics for limited information
#' @param alpha significance level, default alpha=0.05
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @note format for full sample: data[,1] <- study1_treatment; data[,2] <- study1_control; data[,3] <- study2_treatment; data[,4] <- study2_control;...
#' @note format for summary statistics: data[1,] <- c(study1_treatmentlocation, study1_treatmentscale, study1_treatmentsamplesize); data[2,] <- c(study1_controllocation,study1_controlscale,study1_controlsamplesize); data[3,] <- c(study2_treatmentlocation,study2_treatmentscale,study2_treatmentsamplesize);...
#' @note format for summary statistics - two sample hodges-lehmann: data[1,] <- c(study1_shift,study1_scale,study1_treatmentlocation,study1_controllocation,study1_treatmentsamplesize,study1_controlsamplesize);  data[2,] <- c(study2_shift,study2_scale,study2_treatmentlocation,study2_controllocation,study2_treatmentsamplesize,study2_controlsamplesize);...
#' @examples
#' sample1_treat <- rnorm(25,0,1)
#' sample1_con <- rnorm(25,0,1)
#' sample2_treat <- rnorm(25,1,1)
#' sample2_con <- rnorm(25,1,1)
#' dat <- data.frame(sample1_treat,sample1_con,sample2_treat,sample2_con)
#' rmeta(dat, effect = "fixed", method = "MMNP", fullsample = TRUE)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export

rmeta = function(data, effect = "fixed", method = "MMNP", fullsample = TRUE, alpha = 0.05, iternation = 1000){
  a = alpha; n.iter = iternation
  if (effect == "fixed") {
    if (fullsample==TRUE){
      dat_col = ncol(data)
      sumStat = matrix(NA,nrow=dat_col/2,ncol=6)
      sum_MM = matrix(NA,nrow=dat_col/2,ncol=6)
      sum_LS = matrix(NA,nrow=dat_col/2,ncol=6)
      sum_HLKSM = matrix(NA,nrow=dat_col/2,ncol=6)
      sum_WNP = matrix(NA,nrow=dat_col/2,ncol=5)
      sum_WPS = matrix(NA,nrow=dat_col/2,ncol=2)
      sum_WNPLI = matrix(NA,nrow=dat_col/2,ncol=3)
      sum_MMLI = matrix(NA,nrow=dat_col/2,ncol=4)
      sum_HLKSMLI = matrix(NA,nrow=dat_col/2,ncol=4)

      sumStat = getInfo(data,alpha=0.05)
      sum_MM = sumStat[,c(7,8,9,10,1,2)]
      sum_LS = sumStat[,c(3,4,5,6,1,2)]
      sum_HLKSM = sumStat[,c(11,12,13,14,1,2)]
      sum_WNP = sumStat[,c(15,13,14,1,2)]
      sum_WPS = sumStat[,c(15,16)]
      sum_WNPLI = sumStat[,c(15,1,2)]
      sum_MMLI = sumStat[,c(7,8,1,2)]
      sum_HLKSMLI = sumStat[,c(11,12,1,2)]

      if (method == "MMNP") {
        meta_result = FE1S2T(sum_MM,alpha)
      } else if (method == "MMPS") {
        meta_result = FE1SPT(sum_MM,alpha)
      } else if (method == "LSNP") {
        meta_result = FE1S2T(sum_LS,alpha)
      } else if (method == "LSPS") {
        meta_result = FE1SPT(sum_LS,alpha)
      }  else if (method == "HLNP") {
        meta_result = FE1S2T(sum_HLKSM,alpha)
      } else if (method == "HLPS") {
        meta_result = FE1SPT(sum_HLKSM,alpha)
      } else if (method == "WNP") {
        meta_result = FEWNP(sum_WNP,alpha)
      } else if (method == "WPS") {
        meta_result = FEWPS(sum_WPS,alpha)
      }  else if (method == "MMLI") {
        meta_result = FELI(sum_MMLI,alpha)
      } else if (method == "MMLI2") {
        meta_result = FELI2(sum_MMLI,alpha)
      } else if (method == "MMLI3") {
        meta_result = FELI3(sum_MMLI,alpha)
      } else if (method == "HLLI") {
        meta_result = FELI(sum_HLKSMLI,alpha)
      } else if (method == "WNPLI") {
        meta_result = FEWNPLI(sum_WNPLI,alpha)
      }
    } else {
      if (method == "MMNP") {
        meta_result = FE1S2T(data,alpha)
      } else if (method == "MMPS") {
        meta_result = FE1SPT(data,alpha)
      } else if (method == "LSNP") {
        meta_result = FE1S2T(data,alpha)
      } else if (method == "LSPS") {
        meta_result = FE1SPT(data,alpha)
      }  else if (method == "HLNP") {
        meta_result = FE1S2T(data,alpha)
      } else if (method == "HLPS") {
        meta_result = FE1SPT(data,alpha)
      } else if (method == "WNP") {
        meta_result = FEWNP(data,alpha)
      } else if (method == "WPS") {
        meta_result = FEWPS(data,alpha)
      }  else if (method == "MMLI") {
        meta_result = FELI(data,alpha)
      } else if (method == "MMLI2") {
        meta_result = FELI2(data,alpha)
      } else if (method == "MMLI3") {
        meta_result = FELI3(data,alpha)
      } else if (method == "HLLI") {
        meta_result = FELI(data,alpha)
      } else if (method == "WNPLI") {
        meta_result = FEWNPLI(data,alpha)
      }
    }
  }
  else if (effect == "random") {
    if (fullsample==TRUE) {
      dat_col = ncol(data)
      sumStat = matrix(NA,nrow=dat_col/2,ncol=6)
      sum_MM = matrix(NA,nrow=dat_col/2,ncol=6)
      sum_LS = matrix(NA,nrow=dat_col/2,ncol=6)
      sum_HLKSM = matrix(NA,nrow=dat_col/2,ncol=6)
      sum_W = matrix(NA,nrow=dat_col/2,ncol=5)
      #sum_WPS = matrix(NA,nrow=dat_col/2,ncol=2)
      sum_WNPLI = matrix(NA,nrow=dat_col/2,ncol=3)
      sum_MMLI = matrix(NA,nrow=dat_col/2,ncol=4)
      sum_HLKSMLI = matrix(NA,nrow=dat_col/2,ncol=4)

      sumStat = getInfo(data,alpha=0.05)
      sum_MM = sumStat[,c(7,8,9,10,1,2)]
      sum_LS = sumStat[,c(3,4,5,6,1,2)]
      sum_HLKSM = sumStat[,c(11,12,13,14,1,2)]
      sum_W = sumStat[,c(15,13,14,1,2)]
      #sum_WPS = sumStat[,c(15,16)]
      sum_WNPLI = sumStat[,c(15,1,2)]
      sum_MMLI = sumStat[,c(7,8,1,2)]
      sum_HLKSMLI = sumStat[,c(11,12,1,2)]

      if (method == "MMNPREST") {
        meta_result = REST1S2T(sum_MM,alpha=a)
      } else if (method == "MMPSREST") {
        meta_result = REST1SPT(sum_MM,alpha=a)
      } else if (method == "LSNPREST") {
        meta_result = REST1S2T(sum_LS,alpha=a)
      } else if (method == "LSPSREST") {
        meta_result = REST1SPT(sum_LS,alpha=a)
      } else if (method == "HLNPREST") {
        meta_result = REST1S2T(sum_HLKSM,alpha=a)
      } else if (method == "HLPSREST") {
        meta_result = REST1SPT(sum_HLKSM,alpha=a)
      } else if (method == "WNPREST") {
        meta_result = RESTWNP(sum_W,alpha=a)
      } else if (method == "WPSREST") {
        meta_result = RESTWPS(sum_W,alpha=a)
      } else if (method == "MMREBS") {
        meta_result = REBS1S(sum_MMLI,alpha=a)
      } else if (method == "HLREBS") {
        meta_result = REBS1S(sum_HLKSMLI,alpha=a)
      } else if (method == "WREBS") {
        meta_result = REBS2S(WNPLI,alpha=a)
      }
    } else {
      if (method == "MMNPREST") {
        meta_result = REST1S2T(data,alpha=a)
      } else if (method == "MMPSREST") {
        meta_result = REST1SPT(data,alpha=a)
      } else if (method == "LSNPREST") {
        meta_result = REST1S2T(data,alpha=a)
      } else if (method == "LSPSREST") {
        meta_result = REST1SPT(data,alpha=a)
      } else if (method == "HLNPREST") {
        meta_result = REST1S2T(data,alpha=a)
      } else if (method == "HLPSREST") {
        meta_result = REST1SPT(data,alpha=a)
      } else if (method == "WNPREST") {
        meta_result = RESTWNP(data,alpha=a)
      } else if (method == "WPSREST") {
        meta_result = RESTWPS(data,alpha=a)
      } else if (method == "MMREBS") {
        meta_result = REBS1S(data,iternation=n.iter,alpha=a)
      } else if (method == "HLREBS") {
        meta_result = REBS1S(data,iternation=n.iter,alpha=a)
      } else if (method == "WREBS") {
        meta_result = REBS2S(data,iternation=n.iter,alpha=a)
      }
    }
  }
  cl = 1-alpha
  meta_df = data.frame(meta_result,conf.level=cl)
  return(meta_df)
}



