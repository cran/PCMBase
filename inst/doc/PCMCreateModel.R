## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(abind)
library(ape)
library(PCMBase)

if(!requireNamespace("ggtree")) {
  message("Building the vignette requires ggtree R-package. Trying to install.")
  status.ggtree <- try({
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("ggtree", version = "3.9")
  }, silent = TRUE)
  if(class(status.ggtree) == "try-error") {
    stop(
      "The ggtree installation did not succeed. The vignette cannot be built.")
  }
}

## -----------------------------------------------------------------------------
PCMParentClasses.BM_drift <- function(model) {
  c("GaussianPCM", "PCM")
}

## -----------------------------------------------------------------------------
PCMDescribe.BM_drift <- function(model, ...) {
  "Brownian motion model with drift"
}

## -----------------------------------------------------------------------------
PCMCond.BM_drift <- function(
  tree, model, r = 1, metaI = PCMInfo(NULL, tree, model, verbose = verbose),
  verbose=FALSE) {

  
  Sigma_x <- if(is.Global(model$Sigma_x)){as.matrix(model$Sigma_x)}
	     else{as.matrix(model$Sigma_x[,, r])}
  Sigma <- Sigma_x %*% t(Sigma_x)
  if(!is.null(model$Sigmae_x)) {
    Sigmae_x <- if(is.Global(model$Sigmae_x)){as.matrix(model$Sigmae_x)}
		else{as.matrix(model$Sigmae_x[,,r])}
    Sigmae <- Sigmae_x %*% t(Sigmae_x)
  } else {
    Sigmae <- NULL
  }

  if(!is.null(model$h_drift)) { 
    h_drift <- if(is.Global(model$h_drift)) as.vector(model$h_drift) else model$h_drift[, r]
  }else{
    h_drift <- rep(0,nrow(Sigma_x))
  }

  V <- PCMCondVOU(matrix(0, nrow(Sigma), ncol(Sigma)), Sigma, Sigmae)
  omega <- function(t, edgeIndex, metaI) {
    t*h_drift
  }
  Phi <- function(t, edgeIndex, metaI, e_Ht = NULL) {
    diag(nrow(Sigma))
  }
  list(omega = omega, Phi = Phi, V = V)
}

## -----------------------------------------------------------------------------
PCMDescribeParameters.BM_drift <- function(model, ...) {
  list(
    X0 = "trait values at the root",
    h_drift = "drift vector modifying the expectation",
    Sigma_x = "Upper triangular factor of the unit-time variance rate",
    Sigmae_x = "Upper triangular factor of the non-heritable variance 
		or the variance of the measurement error")
}

## -----------------------------------------------------------------------------
PCMListParameterizations.BM_drift <- function(model, ...) {
  list(
    X0 = list(
      c("VectorParameter", "_Global"),
      c("VectorParameter", "_Fixed", "_Global"),
      c("VectorParameter", "_AllEqual", "_Global"),
      c("VectorParameter", "_Omitted")),
    h_drift = list(     
      c("VectorParameter"),
      c("VectorParameter", "_Fixed"),
      c("VectorParameter", "_AllEqual"),
       c("VectorParameter", "_Omitted")),
       
    Sigma_x = list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal")),

    Sigmae_x = list(
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal"),
      c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal", "_Global"),
      c("MatrixParameter", "_Omitted"))
  )
}

## -----------------------------------------------------------------------------
PCMListDefaultParameterizations.BM_drift <- function(model, ...) {
  list(
    X0 = list(
      c("VectorParameter", "_Global"),
      c("VectorParameter", "_Omitted")
    ),
    h_drift = list(     
      c("VectorParameter")),
       
    Sigma_x = list(
        c("MatrixParameter", "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal"),
        c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal"),
        c("MatrixParameter", "_ScalarDiagonal", "_WithNonNegativeDiagonal")
      ),

    Sigmae_x = list(
      c("MatrixParameter", "_Omitted"))
  )
}

## -----------------------------------------------------------------------------
PCMSpecify.BM_drift <- function(model, ...) {
  spec <- list(
    X0 = structure(0.0, class = c('VectorParameter', '_Global'),
                   description = 'trait values at the root'),
    h_drift = structure(0.0, class = c('VectorParameter'),
                   description = 'drift vector modifying the expectation'),
    Sigma_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal', 
					'_WithNonNegativeDiagonal'),
                        description = 'Cholesky factor of the unit-time variance rate'),
    Sigmae_x = structure(0.0, class = c('MatrixParameter', '_UpperTriangularWithDiagonal',
					'_WithNonNegativeDiagonal'),
                         description = 'Upper triangular factor of the non-heritable variance 
                        		or the variance of the measurement error'))
  attributes(spec) <- attributes(model)
  if(is.null(names(spec))) names(spec) <- c('X0', 'h_drift', 'Sigma_x', 'Sigmae_x')
  if(any(sapply(spec, is.Transformable))) class(spec) <- c(class(spec), '_Transformable')
  spec
}

## -----------------------------------------------------------------------------
X0 <- c(5, 2, 1) ## root state
  
## in regime a traits evolve independently
a.Sigma_x <- rbind(c(1.6, 0.0, 0.0),c(0.0, 2.4, 0.0),c(0.0, 0.0, 2.0))
## no jumps at the end of a branch
a.Sigmae_x <- rbind(c(0.0, 0.0, 0.0),c(0.0, 0.0, 0.0),c(0.0, 0.0, 0.0))
a.h_drift<-c(4, 5, 6)

## in regime b evolution is correlated
b.Sigma_x <- rbind(c(1.6, 0.3, 0.3), c(0.0, 0.3, 0.4),c(0.0, 0.0, 2.0))
## no jumps at the end of a branch
b.Sigmae_x <- rbind(c(0.0, 0.0, 0.0),c(0.0, 0.0, 0.0),c(0.0, 0.0, 0.0))
b.h_drift<-c(1, 2, 3)

Sigma_x <- abind(a.Sigma_x, b.Sigma_x, along=3, new.names=list(x=NULL,y=NULL,regime=c('a','b')))
Sigmae_x <- abind(a.Sigmae_x,b.Sigmae_x,along=3,new.names=list(x=NULL,y=NULL,regime=c('a','b')))
h_drift <- abind(a.h_drift, b.h_drift, along=2, new.names=list(xy=NULL, regime=c('a','b')))

PCMBase_model_BM_drift <- PCM("BM_drift", k = 3, regimes = c("a", "b"),
params = list(X0 = X0,h_drift = h_drift[,,drop=FALSE],
Sigma_x = Sigma_x[,,,drop=FALSE],Sigmae_x = Sigmae_x[,,,drop=FALSE]))

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

palette <- PCMColorPalette(2, c("a", "b"))

# Plot the tree with branches colored according to the regimes.
# The following function call works only if the ggtree package is installed, 
# which is not on CRAN:
PCMTreePlot(tree.ab) + ggtree::geom_nodelab(size = 2)

## -----------------------------------------------------------------------------
mData<-PCMSim(tree.ab, PCMBase_model_BM_drift, X0)[,1:N] ## we only want the tip data
## NOTE that observations from different species are in the columns NOT in the rows as 
## in other software

## -----------------------------------------------------------------------------
log_lik<- PCMLik(mData, tree.ab, PCMBase_model_BM_drift)
print(log_lik[1]) ## we just want to print the log-likelihood without the attributes

## -----------------------------------------------------------------------------
## create an vector of appropriate length to store the vectorized model parameters
v_param <- double(PCMParamCount(PCMBase_model_BM_drift))

# load the current model parameters into param
PCMParamLoadOrStore(PCMBase_model_BM_drift, v_param, offset=0, load=FALSE)

print(v_param)

## now create a likelihood function for the particular model and observed data
likFun <- PCMCreateLikelihood(mData, tree.ab, PCMBase_model_BM_drift)

log_lik_from_likFun<-likFun(v_param)
print(log_lik_from_likFun[1])
print(log_lik_from_likFun[1]==log_lik[1])

# modify slightly the model parameters
v_param_2 <- jitter(v_param)

print(v_param_2)

# set the new parameter vector
PCMBase_model_BM_drift_2<-PCMBase_model_BM_drift
PCMParamLoadOrStore(PCMBase_model_BM_drift_2, v_param_2, offset = 0, load=TRUE)

print(PCMBase_model_BM_drift_2)
log_lik_from_likFun_2<-likFun(v_param_2)
print(log_lik_from_likFun_2[1])


