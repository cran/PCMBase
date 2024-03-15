## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ape)
library(PCMBase)

FLAGSuggestsAvailable <- PCMBase::RequireSuggestedPackages()

options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
# scroll to the right in the following listing to see the full model type names 
# and their corresponding alias:
PCMDefaultModelTypes()

## -----------------------------------------------------------------------------
PCMFindMethod(PCMDefaultModelTypes()["A"], "PCMCond")
# The complex maths is implemented in the function PCMCondVOU. You can see its 
# R-code by typing :
# PCMBase::PCMCondVOU

## -----------------------------------------------------------------------------
modelBM <- PCM(model = "BM", k = 2)

## -----------------------------------------------------------------------------
modelBM

## -----------------------------------------------------------------------------
modelBM.ab <- PCM("BM", k = 2, regimes = c("a", "b"))
modelBM.ab

## -----------------------------------------------------------------------------
modelBM.ab$X0[] <- c(5, 2)

## -----------------------------------------------------------------------------
PCMParamCount(modelBM)
PCMParamCount(modelBM.ab)

## -----------------------------------------------------------------------------
# in regime 'a' the traits evolve according to two independent BM processes (starting from the global vecto X0).
modelBM.ab$Sigma_x[,, "a"] <- rbind(c(1.6, 0),
                                  c(0, 2.4))
modelBM.ab$Sigmae_x[,, "a"] <- rbind(c(.1, 0),
                                   c(0, .4))
# in regime 'b' there is a correlation between the traits
modelBM.ab$Sigma_x[,, "b"] <- rbind(c(1.6, .8),
                                  c(.8, 2.4))
modelBM.ab$Sigmae_x[,, "b"] <- rbind(c(.1, 0),
                                   c(0, .4))

## -----------------------------------------------------------------------------
param <- double(PCMParamCount(modelBM.ab))

# load the current model parameters into param
PCMParamLoadOrStore(modelBM.ab, param, offset=0, load=FALSE)

print(param)

# modify slightly the model parameters
param2 <- jitter(param)

print(param2)

# set the new parameter vector
PCMParamLoadOrStore(modelBM.ab, param2, offset = 0, load=TRUE)

print(modelBM.ab)

## ---- results='asis'----------------------------------------------------------
options(digits = 2)
print(PCMTable(modelBM.ab), xtable = TRUE, type="html")

## -----------------------------------------------------------------------------
model.OU.BM <- MixedGaussian(
  k = 3, 
  modelTypes = c(
    BM = paste0("BM",
    "__Omitted_X0",
    "__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x",
    "__Omitted_Sigmae_x"),
    OU = paste0("OU",
    "__Omitted_X0",
    "__H",
    "__Theta",
    "__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x",
    "__Omitted_Sigmae_x")), 
  mapping = c(a = 2, b = 1), 
  Sigmae_x = structure(
    0, 
    class = c("MatrixParameter", "_Omitted"), 
    description = "upper triangular factor of the non-phylogenetic variance-covariance"))

## -----------------------------------------------------------------------------
model.OU.BM$X0[] <- c(NA, NA, NA)
model.OU.BM$`a`$H[,,1] <- cbind(
  c(.1, -.7, .6), 
  c(1.3, 2.2, -1.4), 
  c(0.8, 0.2, 0.9))
model.OU.BM$`a`$Theta[] <- c(1.3, -.5, .2)
model.OU.BM$`a`$Sigma_x[,,1] <- cbind(
  c(1, 0, 0), 
  c(1.0, 0.5, 0), 
  c(0.3, -.8, 1))

model.OU.BM$`b`$Sigma_x[,,1] <- cbind(
  c(0.8, 0, 0), 
  c(1, 0.3, 0), 
  c(0.4, 0.5, 0.3))

## ---- results='asis'----------------------------------------------------------
print(PCMTable(model.OU.BM), xtable = TRUE, type="html")

## -----------------------------------------------------------------------------
# make results reproducible
set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion")

# number of regimes
R <- 2

# number of extant tips
N <- 100

tree.a <- PCMTree(rtree(n=N))
PCMTreeSetLabels(tree.a)
PCMTreeSetPartRegimes(tree.a, part.regime = c(`101` = "a"), setPartition = TRUE)

lstDesc <- PCMTreeListDescendants(tree.a)
splitNode <- names(lstDesc)[which(sapply(lstDesc, length) > N/2 & sapply(lstDesc, length) < 2*N/3)][1]

tree.ab <- PCMTreeInsertSingletons(
  tree.a, nodes = as.integer(splitNode), 
  positions = PCMTreeGetBranchLength(tree.a, as.integer(splitNode))/2)
PCMTreeSetPartRegimes(
  tree.ab,
  part.regime = structure(c("a", "b"), names = as.character(c(N+1, splitNode))), 
  setPartition = TRUE)

## ---- eval=FLAGSuggestsAvailable----------------------------------------------
if(requireNamespace("ggtree")) {
  palette <- PCMColorPalette(2, c("a", "b"))
  
  # Plot the tree with branches colored according to the regimes.
  # The following code works only if the ggtree package is installed, which is not on CRAN. 
  # The tree would not be depicted correctly if ggtree is not installed.
  plTree <- PCMTreePlot(tree.ab)
  plTree <- plTree + ggtree::geom_nodelab(size = 2) 

  plTree
}

## -----------------------------------------------------------------------------
traits <- PCMSim(tree.ab, modelBM.ab, modelBM.ab$X0)

## -----------------------------------------------------------------------------
PCMLik(traits, tree.ab, modelBM.ab)

## -----------------------------------------------------------------------------
# a function of a numerical parameter vector:
likFun <- PCMCreateLikelihood(traits, tree.ab, modelBM.ab)

likFun(param2)

