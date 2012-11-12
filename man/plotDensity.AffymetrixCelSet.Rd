%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotDensity.AffymetrixCelSet}
\alias{plotDensity.AffymetrixCelSet}
\alias{AffymetrixCelSet.plotDensity}
\alias{plotDensity.AffymetrixCelSet}
\alias{plotDensity,AffymetrixCelSet-method}


\title{Plots the densities of all samples}

\description{
 Plots the densities of all samples.
}

\usage{\method{plotDensity}{AffymetrixCelSet}(this, subset=NULL, types=NULL, ..., col=seq_along(this), lty=NULL, lwd=NULL, annotate=TRUE, add=FALSE, verbose=FALSE)}

\arguments{
  \item{subset}{The subset of probes to considered \emph{before} any
    filtering by probe type is applied.
    If a \code{\link[base]{vector}} of \code{\link[base]{double}}s, the cell indices.
    If a scalar \code{\link[base]{double}} in [0,1], the fraction of cells, which can
    be used to speed up the plotting if approximate densities are
    acceptable.
    if \code{\link[base]{NULL}}, all cells are considered.
  }
  \item{types}{The type of probes to include, e.g. \code{"all"},
    \code{"pmmm"}, \code{"pm"}, and \code{"mm"}.}
  \item{...}{Additional arguments passed to
     \code{\link{plotDensity.AffymetrixCelFile}}().}
  \item{col}{A \code{\link[base]{vector}} of colors for each of the arrays.}
  \item{lty}{A \code{\link[base]{vector}} of line types for each of the arrays.}
  \item{lwd}{A \code{\link[base]{vector}} of line widths for each of the arrays.}
  \item{add}{If \code{\link[base:logical]{FALSE}}, a new plot is created, otherwise the generated
    graphics is added to the current plot.}
}

\value{
 Returns nothing.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
