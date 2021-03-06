% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GaussianPCM.R
\name{PCMPLambdaP_1}
\alias{PCMPLambdaP_1}
\title{Eigen-decomposition of a matrix H}
\usage{
PCMPLambdaP_1(H)
}
\arguments{
\item{H}{a numeric matrix}
}
\value{
a list with elements as follows:
\item{lambda}{a vector of the eigenvalues of H}
\item{P}{a squared matrix with column vectors, the eigenvectors of H corresponding to the
eigenvalues in \code{lambda}}
\item{P_1}{the inverse matrix of \code{P}}.
}
\description{
Eigen-decomposition of a matrix H
}
\details{
The function fails with an error message if H is defective, that is, if its matrix of
eigenvectors is computationally singular. The test for singularity is based on the \code{\link{rcond}} function.
}
