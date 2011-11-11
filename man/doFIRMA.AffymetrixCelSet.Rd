%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doFIRMA.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{doFIRMA.AffymetrixCelSet}
\alias{doFIRMA.AffymetrixCelSet}
\alias{AffymetrixCelSet.doFIRMA}
\alias{doFIRMA.AffymetrixCelSet}
\alias{doFIRMA,AffymetrixCelSet-method}

 \alias{doFIRMA.character}
 \alias{doFIRMA}

 \title{Finding Isoforms using Robust Multichip Analysis (FIRMA)}

 \description{
  Finding Isoforms using Robust Multichip Analysis (FIRMA) based on [1].
 }

 \usage{\method{doFIRMA}{AffymetrixCelSet}(csR, ..., flavor=c("v1b", "v1a"), drop=TRUE, verbose=FALSE)}

 \arguments{
  \item{csR}{...}
  \item{...}{Additional arguments passed to \code{\link{FirmaModel}}.}
  \item{flavor}{A \code{\link[base]{character}} string specifying the flavor of FIRMA to use.}
  \item{drop}{If \code{\link[base:logical]{TRUE}}, the FIRMA scores are returned, otherwise
   a named \code{\link[base]{list}} of all intermediate and final results.}
  \item{verbose}{See \code{\link{Verbose}}.}
 }

 \value{
   Returns a named \code{\link[base]{list}}, iff \code{drop == FALSE}, otherwise
   only \code{\link{FirmaSet}} object (containing the FIRMA scores).
 }

 \section{Using a custom exon-by-transcript CDF}{
   It is strongly recommended to use a custom CDF, e.g. "core",
   "extended" or "full" [1].  To use a custom CDF, set it before
   calling this method, i.e. \code{setCdf(csR, cdf)}.
   Do not set the standard "non-supported" Affymetrix CDF
   (see also Section 'Flavors').
 }

 \section{Flavors}{
   If \code{flavor == "v1b"} (default), then the standard
   "non-supported" Affymetrix CDF is used for background correction
   and the quantile normalization steps, and the custom CDF
   is used for the probe summarization and everything that follows.
   The advantage of this flavor is that those two first preprocessing
   steps will remain the same if one later changes to a different
   custom CDF.

   If \code{flavor == "v1a"}, then the custom CDF is used throughout
   all steps of FIRMA, which means that if one changes the custom CDF
   all steps will be redone.
 }

 \references{
  [1] E. Purdom, K. Simpson, M. Robinson, J. Conboy, A. Lapuk & T.P. Speed,
      \emph{FIRMA: a method for detection of alternative splicing from
            exon array data}, Bioinformatics, 2008.\cr
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
