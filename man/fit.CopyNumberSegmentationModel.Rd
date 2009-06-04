%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CopyNumberSegmentationModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fit.CopyNumberSegmentationModel}
\alias{fit.CopyNumberSegmentationModel}
\alias{CopyNumberSegmentationModel.fit}
\alias{fit.CopyNumberSegmentationModel}
\alias{fit,CopyNumberSegmentationModel-method}


\title{Fits the model}

\description{
 Fits the model.
}

\usage{\method{fit}{CopyNumberSegmentationModel}(this, arrays=NULL, chromosomes=getChromosomes(this), force=FALSE, aliased=FALSE, ..., .retResults=FALSE, verbose=FALSE)}

\arguments{
  \item{arrays}{A \code{\link[base]{vector}} of array indices specifying which arrays to
   be considered.  If \code{\link[base]{NULL}}, all are processed.}
  \item{chromosome}{A \code{\link[base]{vector}} of chromosomes indices specifying which
    chromosomes to be considered.  If \code{\link[base]{NULL}}, all are processed.}
  \item{force}{If \code{\link[base:logical]{FALSE}}, the model will not be fitted again if it was
    already fitted.}
  \item{...}{Additional arguments passed to the segmentation method for
    the \code{\link[aroma.core]{RawGenomicSignals}}.}
  \item{.retResults}{If \code{\link[base:logical]{TRUE}}, CBS fit structures are returned for each
    fitted array and chromosome.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a named \code{\link[base]{list}} of named \code{\link[base]{list}}s.
}

\section{Additional arguments to the internal fit function}{
  Arguments in \code{...} are passed down to the internal fit function,
  which means that it is possible to fine tune even further.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{CopyNumberSegmentationModel}}.
}
\keyword{internal}
\keyword{methods}
