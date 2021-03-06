% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCM.R
\name{PCMParentClasses}
\alias{PCMParentClasses}
\title{Parent S3 classes for a model class}
\usage{
PCMParentClasses(model)
}
\arguments{
\item{model}{an S3 object.}
}
\value{
a vector of character string denoting the names of the parent classes
}
\description{
Parent S3 classes for a model class
}
\details{
This S3 generic function is intended to be specified for user models.
This function is called by the `PCM.character` method to determine the parent
classes for a given model class.
}
