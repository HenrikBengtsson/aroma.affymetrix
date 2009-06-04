%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getUnitGroupCellMap.ChipEffectFile}
\alias{getUnitGroupCellMap.ChipEffectFile}
\alias{ChipEffectFile.getUnitGroupCellMap}
\alias{getUnitGroupCellMap.ChipEffectFile}
\alias{getUnitGroupCellMap,ChipEffectFile-method}

\alias{ChipEffectFile.getCellMap}
\alias{getCellMap.ChipEffectFile}
\alias{getCellMap,ChipEffectFile-method}


\title{Gets a (unit, group, cell) index map}

\description{
 Gets a (unit, group, cell) index map.
}

\usage{\method{getUnitGroupCellMap}{ChipEffectFile}(this, units=NULL, force=FALSE, ..., verbose=FALSE)}

\arguments{
  \item{units}{The units for which the map should be returned.
     If \code{\link[base]{NULL}}, all units are considered.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, cached cell indices are ignored.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{data.frame}} with \code{\link[base]{integer}} columns \code{unit}, \code{group},
 and \code{cell}.
}

\examples{\dontrun{
     unit group cell
   1  104     1  335
   2  104     2  336
   3  105     1  337
   4  105     2  338
   5  105     3  339
   6  105     4  340
}}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ChipEffectFile}}.
}
\keyword{internal}
\keyword{methods}
