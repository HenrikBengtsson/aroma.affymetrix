%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{image270.AffymetrixCelFile}
\alias{image270.AffymetrixCelFile}
\alias{AffymetrixCelFile.image270}
\alias{image270.AffymetrixCelFile}
\alias{image270,AffymetrixCelFile-method}


\title{Displays all or a subset of the data spatially}

\description{
 Displays all or a subset of the data spatially.
}

\usage{\method{image270}{AffymetrixCelFile}(this, xrange=c(0, Inf), yrange=c(0, Inf), takeLog=TRUE, interleaved=FALSE, ..., field=c("intensities", "stdvs", "pixels"), col=gray.colors(256), main=getName(this))}

\arguments{
  \item{xrange}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length two giving the left and right
         coordinates of the cells to be returned.}
  \item{yrange}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length two giving the top and bottom
         coordinates of the cells to be returned.}
  \item{...}{Additional arguments passed \code{\link[graphics]{image}} and [...].}
  \item{field}{The data field to be displayed.}
  \item{col}{The color map to be used.}
  \item{main}{The main title of the plot.}
}

\value{
 Returns the (270-degrees rotated) data \code{\link[base]{matrix}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:updateUnits.AffymetrixCelFile]{*updateUnits}()}.
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
