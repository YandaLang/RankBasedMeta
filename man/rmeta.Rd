% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmeta.R
\name{rmeta}
\alias{rmeta}
\title{Rank-based Meta-analysis}
\usage{
rmeta(data, effect = "fixed", method = "MMNP", fullsample = TRUE)
}
\arguments{
\item{data}{input data}

\item{effect}{fixed effects or random effects}

\item{method}{Meta-analysis methods, default MMNP}

\item{fullsample}{data type. full sample or summary statistics for limited information}

\item{alpha}{significance level, default alpha=0.05}

\item{iteration}{number of iterations for bootstrap, only used for random effect bootstrap methods}
}
\description{
Rank-based Meta-analysis
}
\note{
format for full sample: data[,1] <- study1_treatment; data[,2] <- study1_control; data[,3] <- study2_treatment; data[,4] <- study2_control;...

format for summary statistics: data[1,] <- c(study1_treatmentlocation, study1_treatmentscale, study1_treatmentsamplesize); data[2,] <- c(study1_controllocation,study1_controlscale,study1_controlsamplesize); data[3,] <- c(study2_treatmentlocation,study2_treatmentscale,study2_treatmentsamplesize);...

format for summary statistics - two sample hodges-lehmann: data[1,] <- c(study1_shift,study1_scale,study1_treatmentlocation,study1_controllocation,study1_treatmentsamplesize,study1_controlsamplesize);  data[2,] <- c(study2_shift,study2_scale,study2_treatmentlocation,study2_controllocation,study2_treatmentsamplesize,study2_controlsamplesize);...
}
\examples{
sample1_treat <- rnorm(25,0,1)
sample1_con <- rnorm(25,0,1)
sample2_treat <- rnorm(25,1,1)
sample2_con <- rnorm(25,1,1)
dat <- data.frame(sample1_treat,sample1_con,sample2_treat,sample2_con)
rmeta(dat, effect = "fixed", method = "MMNP", fullsample = TRUE)
}
\author{
Yanda Lang, Joseph McKean, Omer Ozturk
}
\keyword{Meta-analysis}
\keyword{robust}
