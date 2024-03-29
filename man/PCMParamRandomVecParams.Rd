% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCMParam.R
\name{PCMParamRandomVecParams}
\alias{PCMParamRandomVecParams}
\title{Generate a random parameter vector for a model using uniform distribution
between its lower and upper bounds.}
\usage{
PCMParamRandomVecParams(
  o,
  k,
  R,
  n = 1L,
  argsPCMParamLowerLimit = NULL,
  argsPCMParamUpperLimit = NULL
)
}
\arguments{
\item{o}{a PCM model object or a parameter}

\item{k}{integer denoting the number of traits.}

\item{R}{integer denoting the number of regimes.}

\item{n}{an integer specifying the number of random vectors to generate}

\item{argsPCMParamLowerLimit, argsPCMParamUpperLimit}{named lists of
arguments passed to
\code{PCMParamLowerLimit} and \code{PCMParamUpperLimit}.}
}
\value{
a numeric matrix of dimension \code{n} x \code{PCMParamCount(o)}.
}
\description{
Generate a random parameter vector for a model using uniform distribution
between its lower and upper bounds.
}
\seealso{
PCMParamLimits PCMParamGetShortVector
}
