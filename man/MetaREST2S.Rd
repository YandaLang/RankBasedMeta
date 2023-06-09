% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetaREST2S.R
\name{MetaREST2S}
\alias{MetaREST2S}
\title{Meta-analysis for random-effect - two-sample traditional}
\usage{
MetaREST2S(info,alpha=0.05)
}
\arguments{
\item{info}{summary statistics}

\item{alpha}{significance level}
}
\description{
Meta-analysis for random-effect - two-sample traditional
}
\note{
format for summary statistics (location, scale, sample size)
}
\examples{
study1 <- c(-80.0, 1243.9)
study2 <- c(-89.0, 1029.3)
info <- rbind(study1,study2)
MetaREST2S(info,alpha=0.05)
}
\author{
Yanda Lang, Joseph McKean
}
\keyword{Meta-analysis}
\keyword{robust}
