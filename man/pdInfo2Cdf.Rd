%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  pdInfo2Cdf.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{pdInfo2Cdf}
\alias{pdInfo2Cdf}

 \alias{PdInfo2Cdf}

 \title{Generates an Affymetrix CDF file from a pdInfo package and a auxillary CEL file for the same chip type}

 \description{
   Generates an Affymetrix CDF file from a pdInfo package and a auxillary CEL file for the same chip type.

   \emph{Disclaimer: This is a user-contributed function.}
 }

 \usage{pdInfo2Cdf(pdpkg, celfile, overwrite=FALSE, verbose=TRUE, ...)}

 \arguments{
  \item{pdpkg}{A \code{\link[base]{character}} string for an existing pdInfo package.}
  \item{celfile}{The pathname to an auxillary CEL for the same chip type.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, an existing CDF is overwritten, otherwise
    an exception is thrown.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
  \item{...}{Not used.}
 }

 \value{
   Returns (invisibly) the pathname to CDF written.
 }

 \author{Maintained by Mark Robinson. Original code by Samuel Wuest.}


\keyword{internal}
