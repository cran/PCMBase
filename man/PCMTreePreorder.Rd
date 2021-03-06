% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCMTree.R
\name{PCMTreePreorder}
\alias{PCMTreePreorder}
\title{Pre-order tree traversal}
\usage{
PCMTreePreorder(tree)
}
\arguments{
\item{tree}{a phylo object with possible singleton nodes (i.e. internal
nodes with one daughter node)}
}
\value{
a vector of indices of edges in tree$edge in pre-order.
}
\description{
Pre-order tree traversal
}
