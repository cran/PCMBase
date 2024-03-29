% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCMTree.R
\name{PCMTreeGetTipsInPart}
\alias{PCMTreeGetTipsInPart}
\title{Get the tips belonging to a part in a tree}
\usage{
PCMTreeGetTipsInPart(tree, part)
}
\arguments{
\item{tree}{a phylo object with an edge.regime member or a PCMTree object}

\item{part}{a character or integer denoting a part name in the tree.}
}
\value{
an integer vector with the ids of the tips belonging to part
}
\description{
Get the tips belonging to a part in a tree
}
\examples{
set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
tree <- ape::rtree(10)
regimes <- sample(letters[1:3], nrow(tree$edge), replace = TRUE)
PCMTreeSetRegimesForEdges(tree, regimes)

if(requireNamespace("ggtree"))
PCMTreePlot(tree) + ggtree::geom_nodelab() + ggtree::geom_tiplab()

part <- PCMTreeGetPartNames(tree)[1]
PCMTreeGetTipsInPart(tree, part)
print(part)

}
\seealso{
\link{PCMTreeGetTipsInRegime}, \link{PCMTreeGetPartNames}, \link{PCMRegimes}, \link{PCMTreeGetPartRegimes}, \link{PCMTreeSetPartRegimes}
}
