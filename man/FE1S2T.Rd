% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FE1S2T.R
\name{FE1S2T}
\alias{FE1S2T}
\title{Meta-analysis for fixed effects non-pooled variance}
\usage{
FE1S2T(info,alpha=0.05)
}
\arguments{
\item{info}{summary statistics}

\item{alpha}{significance level, default alpha=0.05}
}
\description{
Meta-analysis for fixed effects non-pooled variance
}
\note{
format for summary statistics (point estimate 1, point estimate 2, scale 1, scale 2, sample size 1, sample size 2)
}
\examples{
study1 <- c(247.5, 178.0, 13.8, 19.5, 16, 17)
study2 <- c(270.0, 167.0, 9.8, 12.6, 23, 23)
info <- rbind(study1,study2)
FE1S2T(info,alpha=0.05)
}
\author{
Yanda Lang, Joseph McKean, Omer Ozturk
}
\keyword{Meta-analysis}
\keyword{robust}
