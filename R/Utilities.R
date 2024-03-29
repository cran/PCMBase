# Copyright 2016-2019 Venelin Mitov
#
# This file is part of PCMBase.
#
# PCMBase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMBase is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMBase.  If not, see <http://www.gnu.org/licenses/>.


#' True positive rate of a set of binary predictions against their trues
#'
#' @description Let the set of predictions be described by a logical vector
#' `pred`, and let the corresponding trues by described in a logical vector
#' `true` of the same length. Then, the true positive rate is given by the
#' expression:
#' \code{sum(pred & true)/sum(true)}. The false positive rate is given by the
#' expression:
#' \code{sum(pred & !true)/sum(!true)}. If these expressions do not give a
#' finite number, \code{NA_real_} is returned.
#'
#' @param pred,true vectors of the same positive length that can be converted to
#' logical.
#'
#' @return a double between 0 and 1 or \code{NA_real_} if the result is not a
#' finite number.
#' @examples
#' TruePositiveRate(c(1,0,1,1), c(1,1,0,1))
#' TruePositiveRate(c(0,0,0,0), c(1,1,0,1))
#' TruePositiveRate(c(1,1,1,1), c(1,1,0,1))
#' FalsePositiveRate(c(1,0,1,1), c(1,1,0,1))
#' FalsePositiveRate(c(0,0,0,0), c(1,1,0,1))
#' FalsePositiveRate(c(1,1,1,1), c(1,1,0,1))
#' TruePositiveRate(c(1,0,1,1), c(0,0,0,0))
#' FalsePositiveRate(c(1,0,1,1), c(1,1,1,1))
#' @export
TruePositiveRate <- function(pred, true) {
  if(! (length(pred) == length(true)) ) {
    stop("TruePositiveRate:: arguments pred and true must be the same length.")
  }
  pred <- as.logical(pred)
  true <- as.logical(true)
  res <- sum(pred & true)/sum(true)
  if(is.finite(res)) {
    res
  } else {
    NA_real_
  }
}

#' @title True positive rate of a set of binary predictions against their trues
#' @rdname TruePositiveRate
#' @export
FalsePositiveRate <- function(pred, true) {
  pred <- as.logical(pred)
  true <- as.logical(true)
  res <- sum(pred & !true)/sum(!true)
  if(is.finite(res)) {
    res
  } else {
    as.double(NA)
  }
}


#' Check if all packages listed in Suggests are available
#' @return logical TRUE if suggested packages are installed and can be loaded; FALSE
#' otherwise
#' @importFrom utils install.packages
#' @export
RequireSuggestedPackages <- function() {
  isAvailable <- c(
    testthat = requireNamespace("testthat", quietly = TRUE),
    knitr = requireNamespace("knitr", quietly = TRUE),
    rmarkdown = requireNamespace("rmarkdown", quietly = TRUE),
    ggtree  = requireNamespace("ggtree", quietly = TRUE),
    cowplot = requireNamespace("cowplot", quietly = TRUE),
    covr = requireNamespace("covr", quietly = TRUE),
    #mvSLOUCH not called in vignettes
    #mvSLOUCH = requireNamespace("mvSLOUCH", quietly = TRUE),
    BiocManager = requireNamespace("BiocManager", quietly = TRUE)
    )

  if(!sum(isAvailable) == length(isAvailable)) {
    message("PCMBase: The following suggested packages could not be loaded: ",
            toString(names(isAvailable)[!isAvailable]), ". ",
            "The vignette generation and unit tests may fail to execute.")
  }
  sum(isAvailable) == length(isAvailable)
}

#' Check if the PCMBase version corresponds to a dev release
#' @importFrom utils packageDescription
#' @return a logical
#' @export
PCMBaseIsADevRelease <- function() {
  TRUE
}

#' Beautiful model description based on plotmath
#' @description This is an S3 generic that produces a plotmath expression for
#' its argument.
#' @param o a PCM or a parameter object.
#' @param roundDigits an integer, default: 2.
#' @param transformSigma_x a logical indicating if Cholesky transformation should be
#' applied to Cholesky-factor parameters prior to generating the plotmath expression.
#' @return a character string.
#' @export
PCMPlotMath <- function(o, roundDigits = 2, transformSigma_x = FALSE) {
  UseMethod("PCMPlotMath", o)
}

#' @export
PCMPlotMath.default <- function(o, roundDigits = 2, transformSigma_x = FALSE) {
  if(is.numeric(o)) {
    format(round(o, roundDigits), nsmall = roundDigits)
  } else {
    toString(o)
  }
}

#' @export
PCMPlotMath.VectorParameter <- function(o, roundDigits = 2, transformSigma_x = FALSE) {
  if(is.Global(o)) {
    k <- length(o)
    R <- 1
  } else {
    k <- dim(o)[1]
    R <- dim(o)[2]
  }

  if(k <= 2) {
    if(R > 1) {
      res <- "list("
    } else {
      res <- ""
    }
    for(r in seq_len(R)) {
      if(R > 1) {
        res <- paste0(res, "", r, "==")
      }
      res <- paste0(res, 'bgroup("[", ')
      if(k > 1) {
        res <- paste0(res, "atop( ")
      }
      for(row in seq_len(k)) {
        if(is.Global(o)) {
          res <- paste0(res, '"', format(round(o[row], roundDigits), nsmall = roundDigits), '"')
        } else {
          res <- paste0(res, '"', format(round(o[row, r], roundDigits), nsmall = roundDigits), '"')
        }

        if(row != k) {
          res <- paste0(res, ", ")
        }
      }
      if(k > 1) {
        res <- paste0(res, ')')
      }
      res <- paste0(res, ', "]")')
      if(r < R && R > 1) {
        res <- paste0(res, ",")
      }
    }
    if(R > 1) {
      res <- paste0(res, ")")
    }
    res
  } else {
    paste0("[", toString(o), "]")
  }
}

#' @export
PCMPlotMath.MatrixParameter <- function(o, roundDigits = 2, transformSigma_x = FALSE) {
  k <- dim(o)[1]
  if(is.Global(o)) {
    R <- 1
    if(transformSigma_x) {
      if(getOption("PCMBase.Transpose.Sigma_x", FALSE)) {
        o <- t(o) %*% o
      } else {
        o <- o %*% t(o)
      }
    }
  } else {
    R <- dim(o)[3]
    if(transformSigma_x) {
      for(r in seq_len(R)) {
        if(getOption("PCMBase.Transpose.Sigma_x", FALSE)) {
          o[,,r] <- t(o[,,r]) %*% o[,,r]
        } else {
          o[,,r] <- o[,,r] %*% t(o[,,r])
        }
      }
    }
  }

  if(k <= 2) {
    if(R > 1) {
      res <- "list("
    } else {
      res <- ""
    }
    for(r in seq_len(R)) {
      if(R > 1) {
        res <- paste0(res, "", r, "==")
      }
      res <- paste0(res, 'bgroup("[", ')
      if(k > 1) {
        res <- paste0(res, "atop( ")
      }
      for(row in 1:k) {
        res <- paste0(res, "paste(")
        for(col in 1:k) {
          if(is.Global(o)) {
            res <- paste0(res, '"', format(round(o[row, col], roundDigits), nsmall = roundDigits), '"')
          } else {
            res <- paste0(res, '"', format(round(o[row, col, r], roundDigits), nsmall = roundDigits), '"')
          }

          if(col != k) {
            res <- paste0(res, ', " ", ')
          }
        }
        res <- paste0(res, ")")
        if(row != k) {
          res <- paste0(res, ", ")
        }
      }
      if(k > 1) {
        res <- paste0(res, ')')
      }
      res <- paste0(res, ', "]")')

      if(r < R && R > 1) {
        res <- paste0(res, ",")
      }
    }
    if(R > 1) {
      res <- paste0(res, ")")
    }
    res
  } else {
    paste0("[", toString(o), "]")
  }
}

#' @export
PCMPlotMath.PCM <- function(o, roundDigits = 2, transformSigma_x = FALSE) {
  res <- 'bgroup("{", list('
  for(i in seq_along(o)) {
    name <- names(o)[i]
    transformSigma_x <- FALSE
    if(name == "Sigma_x") {
      name <- "Sigma"
      transformSigma_x <- TRUE
    }
    res <- paste0(
      res, name, "==", PCMPlotMath(
        o[[i]], roundDigits = roundDigits,
        transformSigma_x = transformSigma_x))
    if(i < length(o)) {
      res <- paste0(res, ", ")
    }
  }
  res <- paste0(res, '), "}")')
  res
}


#' A fixed palette of n colors
#'
#' @param n an integer defining the number of colors in the resulting palette.
#' @param names a character vector of length `n`.
#' @param colors a vector of n values convertible to colors. Default:
#' \code{structure(hcl(
#' h = seq(15, 375, length = n + 1), l = 65, c = 100)[seq_len(n)],
#' names = names)}
#'
#' @return A vector of character strings which can be used as color
#' specifications by R graphics functions.
#'
#' @importFrom grDevices hcl
#'
#' @export
PCMColorPalette <- function(
  n, names,
  colors = structure(hcl(
    h = seq(15, 375, length = n + 1), l = 65, c = 100)[seq_len(n)],
    names = names)) {

  if(length(colors) != n) {
    stop("colors should be of length n.")
  }
  names(colors) <- names
  colors
}


#' Scatter plot of 2-dimensional data
#' @param X a k x N matrix
#' @param tree a phylo object
#' @param labeledTips a vector of tip-numbers to label (NULL by default)
#' @param palette a named vector of colors
#' @param scaleSizeWithTime logical indicating if the size and the transparency of the points
#'   should reflect the distance from the present (points that are farther away in time with
#'   respect to the present moment, i.e. closer to the root of the tree, are displayed smaller
#'   and more transparent.). By default this is set to \code{!is.ultrametric(tree)}.
#'
#' @param sizeLabels passed to \code{geom_text} to specify the size of tip-labels for the trait-points.
#' @param nudgeLabels a numeric vector of two elements (default: c(0,0)), passed as
#' arguments nudge_x and nudge_y of \code{geom_text}.
#' @param sizePoints,alphaPoints numeric parameters passed as arguments size and alpha to \code{geom_point}.
#' Default: sizePoints = 2, alphaPoints = 1.
#' @param numTimeFacets a number or a numeric vector controlling the creation of different facets
#' corresponding to different time intervals when the tree is non-ultrametric. If a single number,
#' it will be interpreted as an integer specifying the number of facets, each facets corresponding to
#' an equal interval of time. If a numeric vector, it will be used to specify the cut-points for
#' each interval. Default: \code{if(is.ultrametric(tree) || scaleSizeWithTime) 1L else 3}.
#' @param nrowTimeFacets,ncolTimeFacets integers specifying how the time facets should
#' be layed out. Default: \code{nrowTimeFacets = 1L, ncolTimeFacets = numTimeFacets}.
#'
#' @return a ggplot object
#' @importFrom data.table data.table is.data.table setkey := setnames
#' @importFrom ggplot2 ggplot geom_point scale_size_continuous scale_alpha_continuous geom_text aes theme_gray theme scale_color_manual facet_wrap vars
#' @importFrom ape is.ultrametric
#' @export
PCMPlotTraitData2D <- function(
  X, tree,
  sizePoints = 2, alphaPoints = 1,
  labeledTips = NULL, sizeLabels = 8, nudgeLabels = c(0.0, 0.0),
  palette = PCMColorPalette(PCMNumRegimes(tree), PCMRegimes(tree)),
  scaleSizeWithTime = !is.ultrametric(tree),
  numTimeFacets = if(is.ultrametric(tree) || scaleSizeWithTime) 1L else 3L,
  nrowTimeFacets = 1L, ncolTimeFacets = numTimeFacets) {

  tree <- PCMTree(tree)

  # Needed to pass the check
  id <- .N <- time <- regime <- label <- x <- y <- timeFacet <- NULL

  N <- PCMTreeNumTips(tree)
  R <- PCMNumRegimes(tree)

  times <- PCMTreeNodeTimes(tree, tipsOnly = TRUE)
  Xt <- t(X)
  data <- as.data.table(Xt)
  setnames(data, c("x", "y"))
  data[, id:=seq_len(.N)]
  data[, time:=times]
  data[, timeFacet:=if(length(numTimeFacets) > 1L || numTimeFacets > 1L) {
    cut(time, numTimeFacets, include.lowest = TRUE)
  } else {
    factor(paste0("[", toString(range(time)), "]"))
  }]

  data[, regime := factor(PCMTreeGetPartRegimes(tree)[PCMTreeGetPartsForNodes(tree, id)])]
  setkey(data, id)

  if(!is.null(labeledTips)) {
    data[list(labeledTips), label:=id]
  }

  pl <- ggplot(data)
  if(scaleSizeWithTime) {
    # non-ultrametric tree
    pl <- pl + geom_point(aes(x=x, y=y, col=regime, size = time, alpha = time), na.rm = TRUE) +
      scale_size_continuous(range = c(0.2, 2.8)) +
      scale_alpha_continuous(range = c(0.2, 0.75))
  } else {
    pl <- pl + geom_point(
      aes(x=x, y=y, col=regime),
      size = sizePoints, alpha=alphaPoints, na.rm = TRUE)
  }

  pl <- pl + scale_color_manual(name = "regime", values = palette)

  if(!is.null(labeledTips)) {
    pl <- pl + geom_text(
      aes(label = label, x = x, y = y, col = regime),
      size = sizeLabels,
      nudge_x = nudgeLabels[1], nudge_y = nudgeLabels[2],
      show.legend = FALSE,
      check_overlap = TRUE)
  }

  if(length(numTimeFacets) > 1L || numTimeFacets > 1L) {
    pl <- pl + facet_wrap(vars(timeFacet), nrow = nrowTimeFacets, ncol = ncolTimeFacets)
  }

  pl
}

#' A 2D Gaussian distribution density grid in the form of a ggplot object
#' @param mu numerical mean vector of length 2
#' @param Sigma numerical 2 x 2 covariance matrix
#' @param xlim,ylim numerical vectors of length 2
#' @param xNumPoints,yNumPoints integers denoting how many points should the grid contain for each axis.
#' @param ... additional arguments passed to ggplot
#'
#' @return a ggplot object
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom ggplot2 ggplot geom_contour coord_fixed
#'
#' @export
PCMPlotGaussianDensityGrid2D <- function(
  mu, Sigma, xlim, ylim, xNumPoints = 100, yNumPoints = 100, ...) {
  # needed to pass the check
  x <- y <- prob <- NULL

  data.grid <- expand.grid(
    x = seq(xlim[1], xlim[2], length.out = xNumPoints),
    y = seq(ylim[1], ylim[2], length.out = yNumPoints))

  q.samp <- cbind(data.grid, prob = dmvnorm(data.grid, mean = mu, sigma = Sigma))
  ggplot(q.samp, aes(x=x, y=y, z=prob), ...) + coord_fixed(xlim = xlim, ylim = ylim, ratio = 1)
}

#' A 2D sample from Gaussian distribution
#' @param mu numerical mean vector of length 2
#' @param Sigma numerical 2 x 2 covariance matrix
#' @param numPoints an integer denoting how many points should be randomly sampled (see details).
#' @param ... additional arguments passed to ggplot.
#'
#' @return a ggplot object
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom ggplot2 ggplot
#'
#' @details This function generates a random sample of numPoints 2d points using
#' the function rmvnorm from the mvtnorm R-package. Then
#' it produces a ggplot on the generated points.
#' @export
PCMPlotGaussianSample2D <- function(mu, Sigma, numPoints = 1000, ...) {
  # Needed to pass the check
  x <- y <- NULL

  samp <- rmvnorm(numPoints, mean = mu, sigma = Sigma)
  colnames(samp) <- c("x", "y")
  ggplot(as.data.frame(samp), aes(x, y), ...)
}

#' Extract error information from a formatted error message.
#' @param x character string representing the error message.
#' @return Currently the function returns \code{x} unchanged.
#' @export
PCMParseErrorMessage <- function(x) {
  x
}

AsRExpression <- function(v) {
  if(is.list(v)) {
    toString(list(v))
  } else if(is.vector(v)) {
    gsub('^list', 'c', toString(list(as.list(v))), perl=TRUE)
  }
}
