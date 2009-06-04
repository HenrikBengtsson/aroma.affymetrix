%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeLevelModel.fit.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{fit.ProbeLevelModel}
\alias{fit.ProbeLevelModel}
\alias{ProbeLevelModel.fit}
\alias{fit.ProbeLevelModel}
\alias{fit,ProbeLevelModel-method}


 \title{Estimates the model parameters}

 \description{
  Estimates the model parameters for all or a subset of the units.
 }

 \usage{\method{fit}{ProbeLevelModel}(this, units="remaining", ..., force=FALSE, ram=NULL, verbose=FALSE)}

 \arguments{
   \item{units}{The units to be fitted.
     If \code{\link[base]{NULL}}, all units are considered.
     If \code{remaining}, only non-fitted units are considered.
   }
   \item{...}{Arguments passed to \code{\link[aroma.affymetrix:readUnits.ProbeLevelModel]{*readUnits}()}.}
   \item{force}{If \code{\link[base:logical]{TRUE}}, already fitted units are re-fitted, and
     cached data is re-read.}
   \item{ram}{A \code{\link[base]{double}} indicating if more or less units should
     be loaded into memory at the same time.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
  Returns an \code{\link[base]{integer}} \code{\link[base]{vector}} of indices of the units fitted,
  or \code{\link[base]{NULL}} if no units was (had to be) fitted.
 }

 \details{
  All estimates are stored to file.

  The non-array specific parameter estimates together with standard deviation
  estimates and convergence information are stored in one file.

  The parameter estimates specific to each array, typically "chip effects",
  are stored in array specific files.

   Data set specific estimates [L = number of probes]:
    phi [L doubles] (probe affinities), sd(phi) [L doubles],
    isOutlier(phi) [L logicals]

   Algorithm-specific results:
    iter [1 integer], convergence1 [1 logical], convergence2 [1 logical]
    dTheta [1 double]
    sd(eps) - [1 double] estimated standard deviation of the error term

   Array-specific estimates [K = nbr of arrays]:
    theta [K doubles] (chip effects), sd(theta) [K doubles],
    isOutlier(theta) [K logicals]
 
   For each array and each unit group, we store:
     1 theta, 1 sd(theta), 1 isOutlier(theta), i.e. (float, float, bit)
   => For each array and each unit (with \eqn{G_j} groups), we store:
     \eqn{G_j} theta, \eqn{G_j} sd(theta), \eqn{G_j} isOutlier(theta),
   i.e. \eqn{G_j}*(float, float, bit).
   For optimal access we store all thetas first, then all sd(theta) and the
   all isOutlier(theta).
   To keep track of the number of groups in each unit, we have to have a
   (unit, ngroups) map.  This can be obtained from getUnitNames() for the
   AffymetrixCdfFile class.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{ProbeLevelModel}}.
 }


\keyword{internal}
\keyword{methods}
\keyword{IO}
