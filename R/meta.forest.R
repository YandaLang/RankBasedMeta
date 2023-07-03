#' Meta-analysis Forest Plot
#' @description Forest Plot over all studies when full sample is provided
#' @usage
#' meta.forest(data, effect = "fixed", method='MMPS', xlab=c(""), pooled=FALSE)
#' @param data input data
#' @param effect fixed effects or random effects
#' @param method Meta-analysis methods, default MMNP
#' @param xlab label for the plot
#' @param pooled pooled variance or non-pooled variance
#' @param fullsample data type. full sample or summary statistics for limited information
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @examples
#' meta.forest(BMI, effect = "fixed", method='HLNP')
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit metafor
#' @export

meta.forest = function(data, effect="fixed", fullsample=TRUE, method='MMPS', xlab=c(""), pooled=FALSE){
  par(mar = c(2.5, 4, 2, 2))
  if (method=='MMPS' || method=='MMNP' || method=='MMNPREST' || method=='MMPSREST') {
    tempres = rmeta(data, effect = effect, method = method, fullsample = fullsample)
    forest(x=tempres[[2]][,1], vi=tempres[[2]][,2], , header=c("Study","Difference in Medians [95% CI]"), slab=rownames(tempres[[2]]), ylim=c(-1,nrow(tempres[[2]])+3), xlab=xlab)
    abline(h=0)
    addpoly(x=tempres[[1]][1,1], ci.lb=tempres[[1]][1,2], ci.ub=tempres[[1]][1,3], rows=-1,mlab=method)
  } else if (method=='LSPS' || method=='LSNP' || method=='LSNPREST' || method=='LSPSREST') {
    tempres = rmeta(data, effect = effect, method = method, fullsample = fullsample)
    forest(x=tempres[[2]][,1], vi=tempres[[2]][,2], , header=c("Study","Difference in Means [95% CI]"), slab=rownames(tempres[[2]]), ylim=c(-1,nrow(tempres[[2]])+3), xlab=xlab)
    abline(h=0)
    addpoly(x=tempres[[1]][1,1], ci.lb=tempres[[1]][1,2], ci.ub=tempres[[1]][1,3], rows=-1,mlab=method)
  } else if (method=='HLPS' || method=='HLNP' || method=='HLNPREST' || method=='HLPSREST') {
    tempres = rmeta(data, effect = effect, method = method, fullsample = fullsample)
    forest(x=tempres[[2]][,1], vi=tempres[[2]][,2], , header=c("Study","Difference in HL Estimates [95% CI]"), slab=rownames(tempres[[2]]), ylim=c(-1,nrow(tempres[[2]])+3), xlab=xlab)
    abline(h=0)
    addpoly(x=tempres[[1]][1,1], ci.lb=tempres[[1]][1,2], ci.ub=tempres[[1]][1,3], rows=-1,mlab=method)
  } else if (method=='WPS' || method=='WNP' || method=='WPSREST' || method=='WNPREST') {
    tempres = rmeta(data, effect = effect, method = method, fullsample = fullsample)
    forest(x=tempres[[2]][,1], vi=tempres[[2]][,2], , header=c("Study","HL Estimate of the Shift [95% CI]"), slab=rownames(tempres[[2]]), ylim=c(-1,nrow(tempres[[2]])+3), xlab=xlab)
    abline(h=0)
    addpoly(x=tempres[[1]][1,1], ci.lb=tempres[[1]][1,2], ci.ub=tempres[[1]][1,3], rows=-1,mlab=method)
  } else if (method=='MMLI' || method=='MMLI2' || method=='MMLI3' || method=='MMNPREBS' ) {
    if (pooled == TRUE ) {
      tempres = rmeta(data, effect = 'fixed', method = 'MMPS', fullsample = fullsample)
      forest(x=tempres[[2]][,1], vi=tempres[[2]][,2], , header=c("Study","Difference in Medians [95% CI]"), slab=rownames(tempres[[2]]), ylim=c(-1,nrow(tempres[[2]])+3), xlab=xlab)
      abline(h=0)
      tempcombres = rmeta(data, effect = effect, method = method, fullsample = fullsample)
      addpoly(x=tempcombres[[1]][1,1], ci.lb=tempcombres[[1]][1,2], ci.ub=tempcombres[[1]][1,3], rows=-1,mlab=method)
    } else {
      tempres = rmeta(data, effect = 'fixed', method = 'MMNP', fullsample = fullsample)
      forest(x=tempres[[2]][,1], vi=tempres[[2]][,2], , header=c("Study","Difference in Medians [95% CI]"), slab=rownames(tempres[[2]]), ylim=c(-1,nrow(tempres[[2]])+3), xlab=xlab)
      abline(h=0)
      tempcombres = rmeta(data, effect = effect, method = method, fullsample = fullsample)
      addpoly(x=tempcombres[[1]][1,1], ci.lb=tempcombres[[1]][1,2], ci.ub=tempcombres[[1]][1,3], rows=-1,mlab=method)
    }
  } else if (method=='HLLI' || method=='HLNPREBS' ) {
    if (pooled == TRUE ) {
      tempres = rmeta(data, effect = 'fixed', method = 'HLPS', fullsample = fullsample)
      forest(x=tempres[[2]][,1], vi=tempres[[2]][,2], , header=c("Study","Difference in HL Estimates [95% CI]"), slab=rownames(tempres[[2]]), ylim=c(-1,nrow(tempres[[2]])+3), xlab=xlab)
      abline(h=0)
      tempcombres = rmeta(data, effect = effect, method = method, fullsample = fullsample)
      addpoly(x=tempcombres[[1]][1,1], ci.lb=tempcombres[[1]][1,2], ci.ub=tempcombres[[1]][1,3], rows=-1,mlab=method)
    } else {
      tempres = rmeta(data, effect = 'fixed', method = 'HLNP', fullsample = fullsample)
      forest(x=tempres[[2]][,1], vi=tempres[[2]][,2], , header=c("Study","Difference in HL Estimates [95% CI]"), slab=rownames(tempres[[2]]), ylim=c(-1,nrow(tempres[[2]])+3), xlab=xlab)
      abline(h=0)
      tempcombres = rmeta(data, effect = effect, method = method, fullsample = fullsample)
      addpoly(x=tempcombres[[1]][1,1], ci.lb=tempcombres[[1]][1,2], ci.ub=tempcombres[[1]][1,3], rows=-1,mlab=method)
    }
  } else if (method=='WNPLI' || method=='WNPREBS') {
    if (pooled == TRUE ) {
      tempres = rmeta(data, effect = 'fixed', method = 'WNP', fullsample = fullsample)
      forest(x=tempres[[2]][,1], vi=tempres[[2]][,2], , header=c("Study","Difference in HL Estimates [95% CI]"), slab=rownames(tempres[[2]]), ylim=c(-1,nrow(tempres[[2]])+3), xlab=xlab)
      abline(h=0)
      tempcombres = rmeta(data, effect = effect, method = method, fullsample = fullsample)
      addpoly(x=tempcombres[[1]][1,1], ci.lb=tempcombres[[1]][1,2], ci.ub=tempcombres[[1]][1,3], rows=-1,mlab=method)
    }
  }
}



