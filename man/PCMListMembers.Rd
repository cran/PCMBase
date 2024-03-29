% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCM.R
\name{PCMListMembers}
\alias{PCMListMembers}
\title{A vector of access-code strings to all members of a named list}
\usage{
PCMListMembers(
  l,
  recursive = TRUE,
  format = c("$", "$'", "$\\"", "$`", "[['", "[[\\"", "[[`")
)
}
\arguments{
\item{l}{a named list object.}

\item{recursive}{logical indicating if list members should be gone through
recursively. TRUE by default.}

\item{format}{a character string indicating the format for accessing a member.
Acceptable values are \code{c("$", "$'", '$"', '$`', "[['", '[["', '[[`')} of which
the first one is taken as default.}
}
\value{
a vector of character strings denoting each named member of the list.
}
\description{
A vector of access-code strings to all members of a named list
}
\examples{
PCMListMembers(PCMBaseTestObjects$model_MixedGaussian_ab)
PCMListMembers(PCMBaseTestObjects$model_MixedGaussian_ab, format = '$`')
PCMListMembers(PCMBaseTestObjects$tree.ab, format = '$`')
}
