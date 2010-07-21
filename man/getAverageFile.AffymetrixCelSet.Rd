%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.getAverage.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{getAverageFile.AffymetrixCelSet}
\alias{getAverageFile.AffymetrixCelSet}
\alias{AffymetrixCelSet.getAverageFile}
\alias{getAverageFile.AffymetrixCelSet}
\alias{getAverageFile,AffymetrixCelSet-method}


\title{Calculates the mean and the standard deviation of the cell signal (intensity, standard deviation etc.) across the CEL set}

\description{
  Calculates the mean and the standard deviation of the cell signal (intensity, standard deviation etc.) across the CEL set.
}

\usage{\method{getAverageFile}{AffymetrixCelSet}(this, name=NULL, prefix="average", indices="remaining", field=c("intensities", "stdvs"), mean=c("median", "mean"), sd=c("mad", "sd"), na.rm=TRUE, g=NULL, h=NULL, ..., cellsPerChunk=moreCells * 10^7/length(this), moreCells=1, force=FALSE, verbose=FALSE)}

\arguments{
 \item{name}{The label of the calculated parameters.
   If \code{\link[base]{NULL}}, a default name format \code{<prefix>-<mean>-<sd>} is used.}
 \item{indices}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying which cells to consider.
   If \code{"remaining"}, only parameters for cells that have not been
   are calculated.
   If \code{\link[base]{NULL}}, all cells are used.}
 \item{mean}{A \code{\link[base]{character}} of a \code{\link[base]{function}} specifying the function used
   to calculate the average.}
 \item{sd}{A \code{\link[base]{character}} of a \code{\link[base]{function}} specifying the function used
   to calculate the standard deviation.}
 \item{na.rm}{If \code{\link[base:logical]{TRUE}}, \code{\link[base]{NA}}s are excluded before, otherwise not.}
 \item{...}{Not used.}
 \item{cellsPerChunk}{A \code{\link[base]{integer}} specifying the total number of cells
   (across arrays) read into memory per chunk.}
 \item{moreCells}{A \code{\link[base]{double}} scalar indicating if more or less cells
   should be used per chunk.}
 \item{force}{If \code{\link[base:logical]{TRUE}}, parameters for cells already calculated are
   recalculated, otherwise not.}
 \item{verbose}{If \code{\link[base:logical]{TRUE}}, progress details are printed, otherwise not.
   May also be a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns an \code{\link{AffymetrixCelSet}} of the same class as the CEL set
  averaged.
}

\details{
  The parameter estimates are stored as a CEL file of the same class as
  the data files in the set.  The CEL file is named \code{<name>.cel}
  and placed in the directory of the set.
  Currently there is no specific data class for this file, but the average
  cell signals are stored as "intensities", the standard deviation of the
  cell signals as "stddevs", and the number of data points used for each
  estimate is stored as "pixels".
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
