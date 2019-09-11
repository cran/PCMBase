## ----setup---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(abind)
library(PCMBase)

## ------------------------------------------------------------------------
OU <- PCM("OU", k = 3, regimes = c("a", "b"))
class(OU$X0)
class(OU$H)
class(OU$Theta)
class(OU$Sigma_x)
class(OU$Sigmae_x)

## ------------------------------------------------------------------------
OU2 <- PCM(
  paste0(
    "OU_", 
    "_Global_X0_",
    "_Global_H_",
    "_Theta_",
    "_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigma_x_",
    "_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Global_Sigmae_x"), 
  k = 3, regimes = c("a", "b"))
class(OU2$X0)
class(OU2$H)
class(OU2$Theta)
class(OU2$Sigma_x)
class(OU2$Sigmae_x)

## ------------------------------------------------------------------------
BMOU <- MixedGaussian(
  k = 3, 
  modelTypes = c(
    "BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"),
  mapping = c(a = 1, b = 2))

class(BMOU$X0)
class(BMOU$Sigmae_x)

# There are no X0 and Sigmae_x parameters for the regimes because they are omitted:
names(BMOU$a)
names(BMOU$b)

## ------------------------------------------------------------------------
BMOU2 <- MixedGaussian(
  k = 3, 
  modelTypes = c(
    "BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"),
  mapping = c(a = 1, b = 2), 
  Sigmae_x = structure(0, 
    class = c("MatrixParameter",
              "_UpperTriangularWithDiagonal", "_WithNonNegativeDiagonal", "_Fixed", "_Global"),
    description =
      "A fixed upper triangular factor of the non-phylogenetic variance-covariance"))

BMOU2$Sigmae_x[] <- rbind(
    c(0.5, 0.02, 0),
    c(0, 0.2, 0.01),
    c(0, 0, 0.1))

class(BMOU2$X0)
class(BMOU2$Sigmae_x)

# There are no X0 and Sigmae_x parameters for the regimes because they are omitted:
names(BMOU2$a)
names(BMOU2$b)

# Notice that the BMOU model created in the previous example has more parameters than BMOU2:
PCMParamCount(BMOU)
PCMParamCount(BMOU2)

## ------------------------------------------------------------------------
M <- structure(
  rbind(c(0.2, 0.5, 1.2),
        c(0, 0.1, 0.02),
        c(0, 0, 1.02)),
  class = c("MatrixParameter", "_CholeskyFactor", "_Transformable", "_Global"))

Mtransf <- PCMApplyTransformation(M)
Mtransf

## ------------------------------------------------------------------------
# Diagonal
Hs1 <- structure(
  rbind(c(0.2, 0, 0),
        c(0, 0, 0),
        c(0, 0, 1.02)),
  class = c("MatrixParameter", "_Schur", "_Transformable", "_Global"))

PCMApplyTransformation(Hs1)

# Symmetric positive definite with eigenvalues 1.02, 0.1 and 0.02
Hs2 <- structure(
  rbind(c(1.02, 0.5, 1.2),
        c(0, 0.1, 0.02),
        c(0, 0, 0.02)),
  class = c("MatrixParameter", "_Schur", "_Transformable", "_Global"))

PCMApplyTransformation(Hs2)
eigen(PCMApplyTransformation(Hs2))$values

# Asymmetric positive definite with eigenvalues 1.02, 0.1 and 0.02
Hs3 <- structure(
  rbind(c(1.02, 0.5, 1.2),
        c(0.2, 0.1, 0.02),
        c(0.8, 0.1, 0.02)),
  class = c("MatrixParameter", "_Schur", "_Transformable", "_Global"))

PCMApplyTransformation(Hs3)
eigen(PCMApplyTransformation(Hs3))$values


## ------------------------------------------------------------------------
modelObject <- PCM("BM", k = 2L, regimes = c("a", "b", "c"))

# let's assign some values to the model parameters:
vec <- seq_len(PCMParamCount(modelObject))
PCMParamLoadOrStore(modelObject, vec, offset = 0, load=TRUE)

str(modelObject)

## ---- results='asis'-----------------------------------------------------
options(digits = 0)
print(
  PCMTable(modelObject, addTransformed = FALSE, removeUntransformed = FALSE), 
  xtable = TRUE, type='html')

## ------------------------------------------------------------------------
PCMModels(parentClass = "BM")

## ------------------------------------------------------------------------
PCMListParameterizations(structure(0.0, class="BM"))$Sigmae_x

## ------------------------------------------------------------------------
# 1. Filter the list of parametrizations to avoid generating too many S3 methods.
# (note that we could do the same type of filtering for the other parameters).
listParameterizationsBM <- PCMListParameterizations(structure(0.0, class="BM"))
listParameterizationsBM$Sigmae_x <- listParameterizationsBM$Sigmae_x[5]
  
# 2. Generate a table of parametrizations for this list:
dtParameterizations <- PCMTableParameterizations(
  structure(0.0, class="BM"), listParameterizations = listParameterizationsBM)

print(dtParameterizations)

# 3. Generate the parametrizations (optionally, we could select a subset of the
# rows in the data.table)
PCMGenerateParameterizations(structure(0.0, class="BM"), 
                             tableParameterizations = dtParameterizations[])

## ------------------------------------------------------------------------
PCMModels("BM")

## ---- results='asis'-----------------------------------------------------
BMModelGlobalSigmae <- paste0(
  "BM" , 
  "__Global_X0", 
  "__Diagonal_WithNonNegativeDiagonal_Sigma_x", 
  "__Diagonal_WithNonNegativeDiagonal_Global_Sigmae_x")

modelObject2 <- PCM(BMModelGlobalSigmae, k = 2, regimes = c("a", "b", "c"))
vec <- seq_len(PCMParamCount(modelObject2))
PCMParamLoadOrStore(modelObject2, vec, offset = 0, load=TRUE)
    
print(
  PCMTable(modelObject2, addTransformed = FALSE, removeUntransformed = FALSE), 
  xtable = TRUE, type='html')                 

