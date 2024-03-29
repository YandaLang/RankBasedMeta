% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.boxplot.R
\name{meta.boxplot}
\alias{meta.boxplot}
\title{Meta-analysis Boxplot}
\usage{
meta.boxplot(data, horizontal = FALSE)
}
\arguments{
\item{data}{input data}

\item{horizontal}{horizontal or vertical boxplots; default vertical boxes}
}
\description{
Comparison boxplots over all studies
}
\examples{
sample1_treat <- rnorm(25,0,1)
sample1_con <- rnorm(25,0,1)
sample2_treat <- rnorm(25,1,1)
sample2_con <- rnorm(25,1,1)
dat <- data.frame(sample1_treat,sample1_con,sample2_treat,sample2_con)
meta.boxplot(dat)
}
\author{
Yanda Lang, Joseph McKean, Omer Ozturk
}
\keyword{Meta-analysis}
\keyword{robust}
