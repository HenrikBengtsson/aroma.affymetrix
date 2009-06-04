%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CopyNumberChromosomalModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{estimateSds.CopyNumberChromosomalModel}
\alias{estimateSds.CopyNumberChromosomalModel}
\alias{CopyNumberChromosomalModel.estimateSds}
\alias{estimateSds.CopyNumberChromosomalModel}
\alias{estimateSds,CopyNumberChromosomalModel-method}


\title{Estimates the standard deviation of the raw copy numbers (log2-ratios) robustly}

\description{
 Estimates the standard deviation of the raw copy numbers (log2-ratios) robustly using a first-order difference variance estimator, which is
 an estimator that is fairly robust for change points.
}

\usage{\method{estimateSds}{CopyNumberChromosomalModel}(this, arrays=seq(length = nbrOfArrays(this)), chromosomes=getChromosomes(this), ..., verbose=FALSE)}

\arguments{
  \item{arrays}{The arrays to be queried.}
  \item{chromosomes}{The chromosomes to be queried.}
  \item{...}{Additional arguments passed to
     \code{\link[aroma.affymetrix:extractRawCopyNumbers.CopyNumberChromosomalModel]{*extractRawCopyNumbers}()}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a CxK \code{\link[base]{double}} \code{\link[base]{matrix}} where C is the number of chromosomes,
 and K is the number of arrays (arrays are always the last dimension).
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{CopyNumberChromosomalModel}}.
}
\keyword{internal}
\keyword{methods}
