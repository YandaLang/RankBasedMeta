#' Rank-based Meta-analysis
#' @description Rank-based Meta-analysis
#' @usage
#' rmeta(data, effect = "fixed", method = "MM", fullsample = TRUE)
#' @param data input data
#' @param effect fixed effects or random effects
#' @param method Meta-analysis methods
#' @param fullsample data type. full sample or summary statistics
#' @author Yanda Lang, Joseph McKean
#' @note format for full sample: data[,1] <- study1_treatment; data[,2] <- study1_control; data[,3] <- study2_treatment; data[,4] <- study2_control;...
#' @note format for summary statistics: data[1,] <- c(study1_treatmentlocation, study1_treatmentscale, study1_treatmentsamplesize); data[2,] <- c(study1_controllocation,study1_controlscale,study1_controlsamplesize); data[3,] <- c(study2_treatmentlocation,study2_treatmentscale,study2_treatmentsamplesize);...
#' @note format for summary statistics - two sample hodges-lehmann: data[1,] <- c(study1_shift,study1_scale,study1_treatmentlocation,study1_controllocation,study1_treatmentsamplesize,study1_controlsamplesize);  data[2,] <- c(study2_shift,study2_scale,study2_treatmentlocation,study2_controllocation,study2_treatmentsamplesize,study2_controlsamplesize);...
#' @examples
#' sample1_treat <- rnorm(25,0,1)
#' sample1_con <- rnorm(25,0,1)
#' sample2_treat <- rnorm(25,1,1)
#' sample2_con <- rnorm(25,1,1)
#' dat <- data.frame(sample1_treat,sample1_con,sample2_treat,sample2_con)
#' rmeta(dat, effect = "fixed", method = "MM", fullsample = TRUE)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit npsmReg2
#' @export

rmeta = function(data, effect = "fixed", method = "MM", fullsample = TRUE){
  if (effect == "fixed") {
    if (fullsample==TRUE){
      dat_col = ncol(data)
      sumStat = matrix(NA,nrow=dat_col,ncol=10)
      sum_MM = matrix(NA,nrow=dat_col,ncol=4)
      sum_LS = matrix(NA,nrow=dat_col,ncol=4)
      sum_HLKSM = matrix(NA,nrow=dat_col,ncol=4)
      sum_HLMPPT = matrix(NA,nrow=dat_col,ncol=4)
      sum_HL = matrix(NA,nrow=dat_col,ncol=4)
      sum_2SHL = matrix(NA,nrow=dat_col/2,ncol=4)
      sum_2SHL12 = matrix(NA,nrow=dat_col/2,ncol=4)
      sum_2SHL12P = matrix(NA,nrow=dat_col/2,ncol=4)
      sum_2SHLKY = matrix(NA,nrow=dat_col/2,ncol=6)

      sumStat = getInfo(data,alpha=0.05)
      sum_MM = sumStat[,c(5,6,2,1)]
      sum_LS = sumStat[,c(3,4,2,1)]
      sum_HLKSM = sumStat[,c(7,8,2,1)]
      sum_HLMPPT = sumStat[,c(9,10,2,1)]
      sum_HL = sumStat[,c(7,10,2,1)]
      sum_2SHL = get2SHL(data,alpha=0.05)
      sum_2SHL12 = get2SHL12(data,alpha=0.05)
      sum_2SHL12P = get2SHL12P(data,alpha=0.05)
      sum_2SHLKY[,1:2] = sum_2SHL[,1:2]
      for (i in 1:nrow(sum_2SHL)) {
        j = 2*i-1
        sum_2SHLKY[i,3] = sum_MM[j,1]
        sum_2SHLKY[i,4] = sum_MM[j+1,1]
      }
      sum_2SHLKY[,5:6] = sum_2SHL[,3:4]

      if (method == "MM") {
        meta_result = MetaFE1S2T(sum_MM,alpha=0.05)
      } else if (method == "MMPS") {
        meta_result = MetaFE1SPT(sum_MM,alpha=0.05)
      } else if (method == "LSWS") {
        meta_result = MetaFE1S2T(sum_LS,alpha=0.05)
      } else if (method == "LSPS") {
        meta_result = MetaFE1SPT(sum_LS,alpha=0.05)
      } else if (method == "HLKSM") {
        meta_result = MetaFE1S2T(sum_HLKSM,alpha=0.05)
      } else if (method == "HLMPPT") {
        meta_result = MetaFE1SPT(sum_HLMPPT,alpha=0.05)
      } else if (method == "HL2T") {
        meta_result = MetaFE1S2T(sum_HL,alpha=0.05)
      } else if (method == "HLPT") {
        meta_result = MetaFE1SPT(sum_HL,alpha=0.05)
      } else if (method == "2SHL") {
        meta_result = MetaFE2SHL(sum_2SHL,alpha=0.05)
      } else if (method == "2SHL12") {
        meta_result = MetaFE2SHL(sum_2SHL12,alpha=0.05)
      } else if (method == "2SHL12P") {
        meta_result = MetaFE2SHL(sum_2SHL12P,alpha=0.05)
      } else if (method == "C1") {
        meta_result = MetaFEC1(sum_MM,alpha=0.05)
      } else if (method == "C2") {
        meta_result = MetaFEC2(sum_MM,alpha=0.05)
      } else if (method == "C3") {
        meta_result = MetaFEC3(sum_MM,alpha=0.05)
      } else if (method == "D1") {
        meta_result = MetaFED1(sum_HLKSM,alpha=0.05)
      } else if (method == "D2") {
        meta_result = MetaFED2(sum_HLKSM,alpha=0.05)
      } else if (method == "D3") {
        meta_result = MetaFED3(sum_HLKSM,alpha=0.05)
      } else if (method == "E1") {
        meta_result = MetaFEE1(sum_2SHL,alpha=0.05)
      } else if (method == "E2") {
        meta_result = MetaFEE2(sum_2SHL,alpha=0.05)
      } else if (method == "E3") {
        meta_result = MetaFEE3(sum_2SHL,alpha=0.05)
      } else if (method == "KYFM") {
        meta_result = MetaFEKY1S(sum_MM,alpha=0.05)
      } else if (method == "KYFHL") {
        meta_result = MetaFEKY1S(sum_HLKSM,alpha=0.05)
      } else if (method == "KYF2SHL") {
        meta_result = MetaFEKY2S(sum_2SHLKY,alpha=0.05)
      }
    } else {
      if (method == "MM") {
        meta_result = MetaFE1S2T(data,alpha=0.05)
      } else if (method == "MMPS") {
        meta_result = MetaFE1SPT(data,alpha=0.05)
      } else if (method == "LSWS") {
        meta_result = MetaFE1S2T(data,alpha=0.05)
      } else if (method == "LSPS") {
        meta_result = MetaFE1SPT(data,alpha=0.05)
      } else if (method == "HLKSM") {
        meta_result = MetaFE1S2T(data,alpha=0.05)
      } else if (method == "HLMPPT") {
        meta_result = MetaFE1SPT(data,alpha=0.05)
      } else if (method == "HL2T") {
        meta_result = MetaFE1S2T(data,alpha=0.05)
      } else if (method == "HLPT") {
        meta_result = MetaFE1SPT(data,alpha=0.05)
      } else if (method == "2SHL") {
        meta_result = MetaFE2SHL(data,alpha=0.05)
      } else if (method == "2SHL12") {
        meta_result = MetaFE2SHL(data,alpha=0.05)
      } else if (method == "2SHL12P") {
        meta_result = MetaFE2SHL(data,alpha=0.05)
      } else if (method == "C1") {
        meta_result = MetaFEC1(data,alpha=0.05)
      } else if (method == "C2") {
        meta_result = MetaFEC2(data,alpha=0.05)
      } else if (method == "C3") {
        meta_result = MetaFEC3(data,alpha=0.05)
      } else if (method == "D1") {
        meta_result = MetaFED1(data,alpha=0.05)
      } else if (method == "D2") {
        meta_result = MetaFED2(data,alpha=0.05)
      } else if (method == "D3") {
        meta_result = MetaFED3(data,alpha=0.05)
      } else if (method == "E1") {
        meta_result = MetaFEE1(data,alpha=0.05)
      } else if (method == "E2") {
        meta_result = MetaFEE2(data,alpha=0.05)
      } else if (method == "E3") {
        meta_result = MetaFEE3(data,alpha=0.05)
      } else if (method == "KYFM") {
        meta_result = MetaFEKY1S(data,alpha=0.05)
      } else if (method == "KYFHL") {
        meta_result = MetaFEKY1S(data,alpha=0.05)
      } else if (method == "KYF2SHL") {
        meta_result = MetaFEKY2S(data,alpha=0.05)
      }
    }
  } else {
    if (fullsample==TRUE) {
      dat_col = ncol(data)
      sumStat = matrix(NA,nrow=dat_col,ncol=10)
      sum_MM = matrix(NA,nrow=dat_col,ncol=4)
      sum_LS = matrix(NA,nrow=dat_col,ncol=4)
      sum_C4 = matrix(NA,nrow=dat_col,ncol=4)
      sum_D4 = matrix(NA,nrow=dat_col,ncol=4)
      sum_2SHL = matrix(NA,nrow=dat_col/2,ncol=4)
      sum_2SHLKY = matrix(NA,nrow=dat_col/2,ncol=6)
      sum_2SHL12P = matrix(NA,nrow=dat_col/2,ncol=4)

      sumStat = getInfo(data,alpha=0.05)
      sum_MM = sumStat[,c(5,6,2,1)]
      sum_LS = sumStat[,c(3,4,2,1)]
      sum_C4 = sumStat[,c(7,4,2,1)]
      sum_D4 = sumStat[,c(5,4,2,1)]
      sum_2SHL = get2SHL(data,alpha=0.05)
      sum_2SHL12P = get2SHL12P(data,alpha=0.05)
      sum_2SHLKY[,1:2] = sum_2SHL[,1:2]
      for (i in 1:nrow(sum_2SHL)) {
        j = 2*i-1
        sum_2SHLKY[i,3] = sum_MM[j,1]
        sum_2SHLKY[i,4] = sum_MM[j+1,1]
      }
      sum_2SHLKY[,5:6] = sum_2SHL[,3:4]

      if (method == "MMREST") {
        meta_result = MetaREST1S(sum_MM,alpha=0.05)
      } else if (method == "MMREKY") {
        meta_result = MetaREKY1S(sum_MM,alpha=0.05)
      } else if (method == "LSWSREST") {
        meta_result = MetaREST1S(sum_LS,alpha=0.05)
      } else if (method == "2SHLREST") {
        meta_result = MetaREST2S(sum_2SHL,alpha=0.05)
      } else if (method == "2SHLREKY") {
        meta_result = MetaREKY2S(sum_2SHLKY,alpha=0.05)
      } else if (method == "2SHL12PREST") {
        meta_result = MetaREST2S(sum_2SHL12P,alpha=0.05)
      } else if (method == "C4REKY") {
        meta_result = MetaREC4(sum_C4,alpha=0.05)
      } else if (method == "D4REKY") {
        meta_result = MetaRED4(sum_D4,alpha=0.05)
      }
    } else {
      if (method == "MMREST") {
        meta_result = MetaREST1S(data,alpha=0.05)
      } else if (method == "MMREKY") {
        meta_result = MetaREKY1S(data,alpha=0.05)
      } else if (method == "LSWSREST") {
        meta_result = MetaREST1S(data,alpha=0.05)
      } else if (method == "2SHLREST") {
        meta_result = MetaREST2S(data,alpha=0.05)
      } else if (method == "2SHLREKY") {
        meta_result = MetaREKY2S(data,alpha=0.05)
      } else if (method == "2SHL12PREST") {
        meta_result = MetaREST2S(data,alpha=0.05)
      } else if (method == "C4REKY") {
        meta_result = MetaREC4(data,alpha=0.05)
      } else if (method == "D4REKY") {
        meta_result = MetaRED4(data,alpha=0.05)
      }
    }
  }
  meta_df = data.frame(meta_result)
  rownames(meta_df) = c('LB','UB')
  return(meta_df)
}
