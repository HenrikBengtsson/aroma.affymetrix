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

 \title{Generates an Affymetrix CDF file from a Platform Design (PD) package and a auxillary CEL file for the same chip type}

 \description{
   Generates an Affymetrix CDF file from a Platform Design (PD) package and a auxillary CEL file for the same chip type.
   Platform Design (PD) packages are also known as "pdInfo" packages.

   \emph{Disclaimer: This is a user-contributed function.}

   \emph{Instead of using this method, we recommend to use
   \code{\link[=writeCdf.AffyGenePDInfo]{writeCdf}()}
   for the \code{AffyGenePDInfo} class.}
 }

 \usage{
pdInfo2Cdf(pdpkg, celfile, overwrite=FALSE, verbose=TRUE, ...)
}

 \arguments{
  \item{pdpkg}{A \code{\link[base]{character}} string for an existing PD package.}
  \item{celfile}{The pathname to an auxillary CEL for the same chip type.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, an existing CDF is overwritten, otherwise
    an exception is thrown.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
  \item{...}{Not used.}
 }

 \value{
   Returns (invisibly) the pathname to CDF written.
   The CDF filename is generated from the name of the PD package.
 }

 \section{Limitations}{
   The information available in the PD package is limited and does
   not contain all information needed to populate a CDF file.
   In order to workaround these limitations, certain CDF entries
   are set to predefined/hardwired values.
   The 'pbase' and 'tbase' entries of the generated CDF file is
   hardwired to "T" and "A", respectively.  Likewise, the 'groupdirection'
   entry is hardwired to "sense".
 }

 \author{
   Maintained by Mark Robinson.
   Original code by Samuel Wuest.
   Code improvements by Henrik Bengtsson.
 }

 \seealso{
   Instead of using this method, we recommend to use
   \code{\link[=writeCdf.AffyGenePDInfo]{writeCdf}()}
   for the \code{AffyGenePDInfo} class.
 }


\keyword{internal}
