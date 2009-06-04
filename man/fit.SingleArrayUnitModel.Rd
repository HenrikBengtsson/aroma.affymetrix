%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SingleArrayUnitModel.fit.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{fit.SingleArrayUnitModel}
\alias{fit.SingleArrayUnitModel}
\alias{SingleArrayUnitModel.fit}
\alias{fit.SingleArrayUnitModel}
\alias{fit,SingleArrayUnitModel-method}


 \title{Estimates the model parameters}

 \description{
  Estimates the model parameters for all or a subset of the units.
 }

 \usage{\method{fit}{SingleArrayUnitModel}(this, arrays=NULL, units="remaining", ..., force=FALSE, verbose=FALSE)}

 \arguments{
   \item{arrays}{The arrays to be fitted.
     If \code{\link[base]{NULL}}, all arrays are considered.
     If \code{remaining}, only non-fitted arrays are considered.
   }
   \item{units}{The units to be fitted.
     If \code{\link[base]{NULL}}, all units are considered.
     If \code{remaining}, only non-fitted units are considered.
   }
   \item{...}{Arguments passed to \code{\link[aroma.affymetrix:readUnits.SingleArrayUnitModel]{*readUnits}()}.}
   \item{force}{If \code{\link[base:logical]{TRUE}}, already fitted units are re-fitted, and
     cached data is re-read.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
  Returns nothing.
 }

 \details{
  All estimates are stored to file.

  The parameter estimates specific to each array,
  typically "chip effects",
  are stored in array specific files.

   Array-specific estimates [K = nbr of arrays]:
    theta [K doubles] (chip effects), sd(theta) [K doubles],
    isOutlier(theta) [K logicals]
 
   For each array and each unit group, we store:
     1 theta, 1 sd(theta), 1 isOutlier(theta), i.e. (float, float, bit)
   => For each array and each unit (with \eqn{G_j} groups), we store:
     \eqn{G_j} theta, \eqn{G_j} sd(theta), \eqn{G_j} isOutlier(theta),
   i.e. \eqn{G_j}*(float, float, bit).
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{SingleArrayUnitModel}}.
 }


\keyword{internal}
\keyword{methods}
\keyword{IO}
