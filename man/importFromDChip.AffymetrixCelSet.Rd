%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.importFromDChip.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{AffymetrixCelSet$importFromDChip}
\alias{AffymetrixCelSet$importFromDChip}
\alias{importFromDChip.AffymetrixCelSet}
\alias{AffymetrixCelSet.importFromDChip}
\alias{importFromDChip,AffymetrixCelSet-method}

\title{Imports dChip-exported CEL files}

\description{
 Imports dChip-exported CEL files into a directory structure recognized by this package.
 ASCII CEL files are converted to binary CEL files, and for chip types
 where the array data is rotated 90-degrees counter clockwise by dChip,
 the data is rotated back.

 As of 2007-03-28, dChip rotates data for exon, tiling, and
 Mapping 500K arrays.
}

\usage{
## Static method (use this):
## AffymetrixCelSet$importFromDChip(path, name=NULL, tags="dChip", rootPath="probeData", rotateBack=NA, ..., skip=TRUE, verbose=FALSE)

## Don't use the below:
\method{importFromDChip}{AffymetrixCelSet}(static, path, name=NULL, tags="dChip", rootPath="probeData", rotateBack=NA, ..., skip=TRUE, verbose=FALSE)
}

\arguments{
  \item{path}{The path to all dChip-exported CEL files.}
  \item{name}{The name of the output data set.
    If \code{\link[base]{NULL}}, the name is inferred from the source path.}
  \item{tags}{Tags added to the imported data set.}
  \item{rootPath}{The root path where to store the data set.}
  \item{rotateBack}{If \code{\link[base:logical]{TRUE}}, the dChip-rotated array data is rotated
    back. If \code{\link[base]{NA}}, this is inferred from the chip type name.}
  \item{...}{Additional arguments passed to \code{byPath()}.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, already converted files are not re-converted.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link{AffymetrixCelSet}}.
}

\details{
 Note that dChip stores CEL intensities as 16-bit integers in its *.dcp
 files, although CEL files are capable or holding floats.
 This means that you might loose precision if first import data in to
 dChip and then export it data again.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
