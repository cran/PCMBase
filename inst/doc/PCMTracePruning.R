## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ape)
library(PCMBase)
library(data.table)
library(xtable)

FLAGSuggestsAvailable <- PCMBase::RequireSuggestedPackages()

options(digits = 2)
# specify either 'html' or 'latex'
tableOutputType <- 'html' # 'latex' is used for generating the S-tables in the ms.

## ----treeAndData--------------------------------------------------------------
library(ape); 
library(PCMBase);

# Non-ultrametric phylogenetic tree of 5 tips in both examples:
treeNewick <- "((5:0.8,4:1.8)7:1.5,(((3:0.8,2:1.6)6:0.7)8:0.6,1:2.6)9:0.9)0;"
tree <- PCMTree(read.tree(text = treeNewick))
# Partitioning the tree in two parts and assign the regimes:
PCMTreeSetPartRegimes(tree, part.regime = c(`6`=2), setPartition = TRUE, inplace = TRUE)

pOrder <- c(PCMTreeGetLabels(tree)[tree$edge[PCMTreePostorder(tree), 2]], "0")

# Trait-data:
X <- cbind(
  c(0.3, NaN, 1.4), 
  c(0.1, NaN, NA), 
  c(0.2, NaN, 1.2), 
  c(NA, 0.2, 0.2), 
  c(NA, 1.2, 0.4))

colnames(X) <- as.character(1:5)

## ----PlotTreeAndData, include=FALSE, eval=FALSE-------------------------------
#  library(tikzDevice); library(ggplot2); library(data.table);
#  # 4. Plotting the tree, the data and the active coordinate vectors:
#  tipValueLabels <- data.table(
#    node = seq_len(PCMTreeNumTips(tree)),
#    valueLabel = paste0(
#      "$\\vec{x}_{", tree$tip.label, "}=(", apply(X[, tree$tip.label], 2, toString), ")^T$"),
#    parse = TRUE)
#  
#  # Determine the active coordinates for X:
#  k_i <- PCMPresentCoordinates(X[, tree$tip.label], tree, NULL)
#  
#  
#  dtNodes <- PCMTreeDtNodes(tree)
#  dtNodes[, kLabel:=paste0(
#    "$\\vec{k}_{", endNodeLab, "}=(", sapply(endNode, function(i) toString(which(k_i[,i]))), ")^T$")]
#  dtNodes[, kLabel2:=kLabel]
#  dtNodes[endNodeLab == "8", kLabel:=NA]
#  dtNodes[endNodeLab != "8", kLabel2:=NA]
#  dtNodes[, tLabel:=paste0("$t_{", endNodeLab, "}=", endTime-startTime, "$")]
#  dtNodes[endNodeLab == "0", tLabel:=NA]
#  
#  tikz(file = "TreeMGPMExample.tex", width = 8, height = 5)
#  palette <- PCMColorPalette(2, names = c("1", "2"), colors = c("black", "orange"))
#  
#  plTree <- PCMTreePlot(tree, palette = palette, size=2)
#  
#  # Plot the tree with branches colored according to the regimes.
#  # The following code works correctly only if the ggtree package is installed,
#  # which is not on CRAN.
#  if(requireNamespace("ggtree")) {
#  
#    plTree <- plTree +
#      ggtree::geom_nodelab(geom = "label", color = "red") +
#      ggtree::geom_tiplab(geom = "label", color = "black")
#  
#    plTree <- plTree %<+% tipValueLabels %<+% dtNodes[, list(node = endNode, kLabel, kLabel2, tLabel)]
#  
#    plTree <- plTree +
#      ggtree::geom_tiplab(geom = "text", aes(label = valueLabel), color = "black", hjust = -0.2, vjust = -1.1) +
#      ggtree::geom_tiplab(geom = "text", aes(label = kLabel), color = "black", hjust = -0.4, vjust = 1.1) +
#      ggtree::geom_nodelab(geom = "text", aes(label = kLabel), color = "red", hjust = -0.2, vjust = 0.6) +
#      ggtree::geom_nodelab(geom = "text", aes(label = kLabel2), color = "red", hjust = 0.4, vjust = 2.8) +
#      geom_text(aes(x = branch, label = tLabel), vjust = -0.8, color = "black") +
#      scale_x_continuous(limits = c(0, 5.2)) + scale_y_continuous(limits = c(0.8, 5.2))
#  }
#  
#  plTree
#  
#  dev.off()

## ---- results='asis'----------------------------------------------------------
model.OU.BM <- MixedGaussian(
  k = nrow(X), 
  modelTypes = c(
    BM = "BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
    OU = "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"), 
  mapping = c(2, 1), 
  Sigmae_x = structure(
    0, 
    class = c("MatrixParameter", "_Omitted", 
              description = "upper triangular factor of the non-phylogenetic variance-covariance")))

model.OU.BM <- PCMApplyTransformation(model.OU.BM)
model.OU.BM$X0[] <- c(NA, NA, NA)
model.OU.BM$`1`$H[,,1] <- cbind(
  c(.1, -.7, .6), 
  c(1.3, 2.2, -1.4), 
  c(0.8, 0.2, 0.9))
model.OU.BM$`1`$Theta[] <- c(1.3, -.5, .2)
model.OU.BM$`1`$Sigma_x[,,1] <- cbind(
  c(1, 0, 0), 
  c(1.0, 0.5, 0), 
  c(0.3, -.8, 1))

model.OU.BM$`2`$Sigma_x[,,1] <- cbind(
  c(0.8, 0, 0), 
  c(1, 0.3, 0), 
  c(0.4, 0.5, 0.3))

print(
  PCMTable(model.OU.BM, removeUntransformed = FALSE), 
  xtable = TRUE, type=tableOutputType)

## ----add-to-PCMBaseTestObjects, include = FALSE, eval=FALSE-------------------
#  # add these objects to the PCMBaseTestObjects (needed for the coding examples).
#  PCMBaseTestObjects[["tree"]] <- tree
#  PCMBaseTestObjects[["X"]] <- X[, tree$tip.label]
#  PCMBaseTestObjects[["model.OU.BM"]] <- model.OU.BM
#  
#  usethis::use_data(PCMBaseTestObjects, overwrite = TRUE)

## ---- echo = TRUE-------------------------------------------------------------
options(digits = 4)
# Variant 1: 
PCMLik(X[, tree$tip.label], tree, model.OU.BM)

# Variant 2: First we call the function PCMInfo to obtain a meta-information object.
metaI.variant2 <- PCMInfo(X[, tree$tip.label], tree, model.OU.BM)
# Then, we manually change the vector of present coordinates for the root node.
# The pc-matrix is a k x M matrix of logical values, each column corresponding
# to a node. The active coordinates are indicated by the TRUE entries.
# To prevent assigning to the wrong column in the pc-table, we first assign
# the node-labels as column nanmes.
colnames(metaI.variant2$pc) <- PCMTreeGetLabels(tree)
metaI.variant2$pc[, "0"] <- c(TRUE, FALSE, TRUE)
# After the change, the pc-matrix looks like this:
metaI.variant2$pc
# And the log-likelihood value is:
PCMLik(X[, tree$tip.label], tree, model.OU.BM, metaI = metaI.variant2)

# Variant 3: We set all NaN values in X to NA, to indicate that these are
# missing measurements
X3 <- X
X3[is.nan(X3)] <- NA_real_
PCMLik(X3[, tree$tip.label], tree, model.OU.BM)

## ----traceTable1R-------------------------------------------------------------
traceTable1 <- PCMLikTrace(X[, tree$tip.label], tree, model.OU.BM)
traceTable1[, node:=.I]
setkey(traceTable1, i)

## ----traceTable2R-------------------------------------------------------------
traceTable2 <- PCMLikTrace(
  X[, tree$tip.label], tree, model.OU.BM, metaI = metaI.variant2)
traceTable2[, node:=.I]
setkey(traceTable2, i)

## ----traceTable3R-------------------------------------------------------------
traceTable3 <- PCMLikTrace(X3[, tree$tip.label], tree, model.OU.BM)
traceTable3[, node:=.I]
setkey(traceTable3, i)

## ----omegaPhiVOU--------------------------------------------------------------
# OU parameters:
H <- model.OU.BM$`1`$H[,,1]
theta <- model.OU.BM$`1`$Theta[,1]
Sigma <- model.OU.BM$`1`$Sigma_x[,,1] %*% t(model.OU.BM$`1`$Sigma_x[,,1])

# Eigenvalues of H: these can be complex numbers
lambda <- eigen(H)$values

# Matrix of eigenvectors of H: again, these can be complex
P <- eigen(H)$vectors
P_1 <- solve(P)

# vectors of active coordinates:
pc <- PCMInfo(X[, tree$tip.label], tree, model.OU.BM)$pc

# length of the branch leading to tip 1 (2.6):
t1 <- PCMTreeDtNodes(tree)[endNodeLab == "1", endTime - startTime]

# active coordinates for tip 1 and its parent:
k1 <- pc[, match("1", PCMTreeGetLabels(tree))]
k9 <- pc[, match("9", PCMTreeGetLabels(tree))]

# k x k matrix formed from the pairs of lambda-values and t1 (see Eq. 19):
LambdaMat <- matrix(0, 3, 3)
for(i in 1:3) 
  for(j in 1:3) 
    LambdaMat[i,j] <- 1/(lambda[i]+lambda[j])*(1-exp(-(lambda[i]+lambda[j])*t1))

# omega, Phi, V for tip 1:
print(omega1 <- (diag(1, 3, 3)[k1, ] - expm::expm(-H*t1)[k1, ]) %*% theta[])
print(Phi1 <- expm::expm(-H*t1)[k1, k9])
print(V1 <- (P %*% (LambdaMat * (P_1 %*% Sigma %*% t(P_1))) %*% t(P))[k1, k1])

## ----omegaPhiVBM--------------------------------------------------------------
# BM parameter:
Sigma <- model.OU.BM$`2`$Sigma_x[,,1] %*% t(model.OU.BM$`2`$Sigma_x[,,1])

# vectors of active coordinates:
pc <- PCMInfo(X[, tree$tip.label], tree, model.OU.BM)$pc

# length of the branch leading to tip 2 (1.6):
t2 <- PCMTreeDtNodes(tree)[endNodeLab == "2", endTime - startTime]

# active coordinates for tip 1 and its parent:
k2 <- pc[, match("2", PCMTreeGetLabels(tree))]
k6 <- pc[, match("6", PCMTreeGetLabels(tree))]

# omega, Phi, V for tip 1:
print(omega2 <- as.matrix(rep(0, 3)[k2]))
print(Phi2 <- as.matrix(diag(1, 3, 3)[k2, k6]))
print(V2 <- as.matrix((t2*Sigma)[k2, k2]))

## ----omegaPhiV1R, results='asis'----------------------------------------------
options(digits = 2)
cat(FormatTableAsLatex(
  traceTable1[list(pOrder), list(j, i, t_i, k_i, omega_i, Phi_i, V_i, V_1_i)], 
  type = tableOutputType))

## ----omegaPhiV2R, results='asis'----------------------------------------------
cat(FormatTableAsLatex(
  traceTable2[list(pOrder), list(j, i, t_i, k_i, omega_i, Phi_i, V_i, V_1_i)], 
  type = tableOutputType))

## ----omegaPhiV3R, results='asis'----------------------------------------------
options(digits = 2)
cat(FormatTableAsLatex(
  traceTable3[list(pOrder), list(j, i, t_i, k_i, omega_i, Phi_i, V_i, V_1_i)], 
  type = tableOutputType))

## ----AbCdEf-------------------------------------------------------------------
# For tip 1. We directly apply Eq. 2, Thm 1:
# We can safely use the real part of V1 (all imaginary parts are 0):
print(V1)
V1 <- Re(V1)
V1_1 <- solve(V1)

print(A1 <- -0.5*V1_1)
print(E1 <- t(Phi1) %*% V1_1)
print(b1 <- V1_1 %*% omega1)
print(C1 <- -0.5 * E1 %*% Phi1)
print(d1 <- -E1 %*% omega1)
print(f1 <- -0.5 * (t(omega1) %*% V1_1 %*% omega1 + sum(k1)*log(2*pi) + log(det(V1))))

## ----AbCdEf1R, results='asis'-------------------------------------------------
cat(FormatTableAsLatex(
  traceTable1[list(pOrder), list(j, i, k_i, A_i, b_i, C_i, d_i, E_i, f_i)], 
  type = tableOutputType))

## ----AbCdEf2R, results='asis'-------------------------------------------------
cat(FormatTableAsLatex(
  traceTable2[list(pOrder), list(j, i, k_i, A_i, b_i, C_i, d_i, E_i, f_i)], 
  type = tableOutputType))

## ----AbCdEf3R, results='asis'-------------------------------------------------
cat(FormatTableAsLatex(
  traceTable3[list(pOrder), list(j, i, k_i, A_i, b_i, C_i, d_i, E_i, f_i)], 
  type = tableOutputType))

## ----LmrVariant1--------------------------------------------------------------
# For tip 2 with parent node 6, we use the following terms stored in Table S5:
A2 <- matrix(-0.17)
b2 <- 0.0
C2 <- rbind(c(-0.17, 0), 
            c(0, 0))
d2 <- c(0.0, 0.0)
E2 <- matrix(c(0.35, 0), nrow = 2, ncol = 1)
f2 <- -1.45
k2 <- 1

# Now we apply Eq. S3:
print(L62 <- C2)
print(m62 <- d2 + E2 %*% X[k2, "2", drop = FALSE])
print(r62 <- t(X[k2, "2", drop = FALSE]) %*% A2 %*% X[k2, "2", drop = FALSE] + 
        t(X[k2, "2", drop = FALSE]) %*% b2 + f2)

# For tip 3 with parent node 6, applying Eq. S3, we obtain (see Table S8):
L63 <- rbind(c(-0.38, 0.51),
             c(0.51, -7.62))
m63 <- c(-1.07, 18.09)
r63 <- -11.41

# Now, we sum the terms L6i, m6i and r6i over all daughters of 6 (i) to obtain:
print(L6 <- L62 + L63)
print(m6 <- m62 + m63)
print(r6 <- r62 + r63)

# Using Eq. S3, we obtain L86, m86, r86, using the values for A,b,C,d,E,f in Table S5:
A6 <- rbind(c(-0.44, 0.58),
            c(0.58, -8.71))
b6 <- c(0.0, 0.0)
C6 <- rbind(c(-0.44, 0.58),
            c(0.58, -8.71))
d6 <- c(0.0, 0.0)
E6 <- rbind(c(0.87, -1.16),
            c(-1.16, 17.42))
f6 <- -0.52
k6 <- c(1, 3)

print(L86 <- C6 - (1/4)*E6 %*% solve(A6 + L6) %*% t(E6))
print(m86 <- d6 - (1/2)*E6 %*% solve(A6 + L6) %*% (b6+m6))
print(r86 <- f6+r6+(length(k6)/2)*log(2*pi) - 
        (1/2)*log(det(-2*(A6+L6))) -
        (1/4)*t(b6+m6) %*% solve(A6+L6) %*% (b6+m6))

# Because 8 is a singleton node, we immediately obtain L8, m8, r8:
L8 <- L86; m8 <- m86; r8 <- r86;

## ----Lmr1R, results='asis'----------------------------------------------------
options(digits = 3)
cat(FormatTableAsLatex(
  traceTable1[
    list(pOrder), 
    list(
      j, i, X_i, k_i, 
      L_i, m_i, r_i,
      `L_{ji}`, `m_{ji}`, `r_{ji}`)], 
  
  type = tableOutputType))

## ----Lmr2R, results='asis'----------------------------------------------------
cat(FormatTableAsLatex(
  traceTable2[
    list(pOrder), 
    list(
      j, i, X_i, k_i, 
      L_i, m_i, r_i,
      `L_{ji}`, `m_{ji}`, `r_{ji}`)], 
  
  type = tableOutputType))

## ----Lmr3R, results='asis'----------------------------------------------------
cat(FormatTableAsLatex(
  traceTable3[
    list(pOrder), 
    list(
      j, i, X_i, k_i, 
      L_i, m_i, r_i, 
      `L_{ji}`, `m_{ji}`, `r_{ji}`)], 
  
  type = tableOutputType))

## -----------------------------------------------------------------------------
# Variant 1.
# Copy the values of L0, m0 and r0 from Table S8:
L0 <- rbind(c(-0.192, 0.214, 0.178),
            c(0.214, -0.313, -0.265),
            c(0.178, -0.265, -0.230))
m0 <- c(0.96, 0.026, 0.255)
r0 <- -18.377

# Use Eq. S2 to estimate the optimal X0:
print(t(x0Hat <- -0.5*solve(L0) %*% m0))
# Use Eq. S1 to calculate the log-likelihood:
print(ll0 <- t(x0Hat) %*% L0 %*% x0Hat + t(x0Hat) %*% m0 + r0)

# Variant 2.
# Copy the values of L0, m0 and r0 from Table S9:
L0 <- rbind(c(-0.192, 0.178),
            c(0.178, -0.230))
m0 <- c(0.96, 0.255)
r0 <- -18.377

# Use Eq. S2 to estimate the optimal X0:
print(t(x0Hat <- -0.5*solve(L0) %*% m0))

# Use Eq. S1 to calculate the log-likelihood:
print(ll0 <- t(x0Hat) %*% L0 %*% x0Hat + t(x0Hat) %*% m0 + r0)

# Variant 3.
# The function PCMLikTrace generates a data.table with the values of 
# omega, Phi, V, A, b, C, d, E, f, L, m, r. 
traceTable3 <- PCMLikTrace(X3[, tree$tip.label], tree, model.OU.BM)
# The column i corresponds to the node label in the tree as depicted on Fig. 1:
setkey(traceTable3, i)

options(digits = 4)
# Variant 3.
# Copy the values of L0, m0 and r0 from the traceTable object (these values have
# the maximal double floating point precision):
print(L0 <- traceTable3[list("0")][["L_i"]][[1]])
print(m0 <- traceTable3[list("0")][["m_i"]][[1]])
print(r0 <- traceTable3[list("0")][["r_i"]][[1]])

# Notice the exact match with the values for variant 3 reported in Fig. S5:
print(t(x0Hat <- -0.5*solve(L0) %*% m0))
print(ll0 <- t(x0Hat) %*% L0 %*% x0Hat + t(x0Hat) %*% m0 + r0)

