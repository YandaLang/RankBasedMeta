% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetaFEE2.R
\name{MetaFEE2}
\alias{MetaFEE2}
\title{Meta-analysis for fixed-effect - E2}
\usage{
MetaFEE2(info,alpha=0.05)
}
\arguments{
\item{info}{summary statistics}

\item{alpha}{significance level}
}
\description{
Meta-analysis for fixed-effect - E2
}
\note{
format for info (2SHL,sample size 1,sample size 2)
}
\examples{
study1 <- c(-80.0, 17, 16)
study2 <- c(-89.0, 23, 23)
info <- rbind(study1,study2)
MetaFEE2(info,alpha=0.05)
}
\author{
Yanda Lang, Joseph McKean
}
\keyword{Meta-analysis}
\keyword{robust}
