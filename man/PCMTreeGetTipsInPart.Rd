% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCMTree.R
\name{PCMTreeGetTipsInPart}
\alias{PCMTreeGetTipsInPart}
\title{Get the tips belonging to a part in a tree}
\usage{
PCMTreeGetTipsInPart(tree, part)
}
\arguments{
\item{tree}{a phylo object with an edge.part member}

\item{part}{a character or integer belonging to tree$edge.part}
}
\value{
an integer vector with the ids of the tips belonging to part
}
\description{
Get the tips belonging to a part in a tree
}