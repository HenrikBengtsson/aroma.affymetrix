%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotDensity.AffymetrixCelFile}
\alias{plotDensity.AffymetrixCelFile}
\alias{AffymetrixCelFile.plotDensity}
\alias{plotDensity,AffymetrixCelFile-method}

\title{Plots the density of the probe signals on the array}

\description{
 Plots the density of the probe signals on the array.
}

\usage{
\method{plotDensity}{AffymetrixCelFile}(this, subset=NULL, types=NULL, ..., xlim=c(0, 16), xlab=NULL, ylab="density (integrates to one)", log=TRUE, annotate=TRUE, verbose=FALSE)
}

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
             \code{\link[aroma.light]{plotDensity.numeric}}.}
  \item{xlim}{The range on the x axis.}
  \item{xlab,ylab}{The labels on the x and the y axes.}
  \item{log}{If \code{\link[base:logical]{TRUE}}, the density of the log (base 2) values are
     used, otherwise the non-logged values.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
}

\value{
 Returns nothing.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
