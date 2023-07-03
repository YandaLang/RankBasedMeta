#' Meta-analysis Boxplot
#' @description Comparison boxplots over all studies
#' @usage
#' meta.boxplot(data, horizontal = FALSE)
#' @param data input data
#' @param horizontal horizontal or vertical boxplots; default vertical boxes
#' @author Yanda Lang, Joseph McKean, Omer Ozturk
#' @examples
#' sample1_treat <- rnorm(25,0,1)
#' sample1_con <- rnorm(25,0,1)
#' sample2_treat <- rnorm(25,1,1)
#' sample2_con <- rnorm(25,1,1)
#' dat <- data.frame(sample1_treat,sample1_con,sample2_treat,sample2_con)
#' meta.boxplot(dat)
#' @keywords robust
#' @keywords Meta-analysis

#' @import stats Rfit
#' @export

meta.boxplot = function(data, horizontal = FALSE){
  par(mar = c(2.5, 4, 2, 2))
  boxplot(data, horizontal = horizontal)
}



