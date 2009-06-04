%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnChipEffectSet.importFromApt.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{importFromApt.CnChipEffectSet}
\alias{importFromApt.CnChipEffectSet}
\alias{CnChipEffectSet.importFromApt}
\alias{importFromApt.CnChipEffectSet}
\alias{importFromApt,CnChipEffectSet-method}


\title{Imports copy-number estimates from an APT summary file}

\description{
 Imports copy-number estimates from an APT summary file.
 Currently only total copy-number estimates can be imported.
}

\synopsis{importFromApt.CnChipEffectSet(static, filename, path=NULL, combineAlleles=TRUE, cdf, ..., skip=TRUE, verbose=FALSE)}
\usage{CnChipEffectSet$importFromApt(filename, path=NULL, combineAlleles=TRUE, cdf, ..., skip=TRUE, verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the APT summary file.}
  \item{path}{An optional path to the file.}
  \item{combineAlleles}{If \code{\link[base:logical]{TRUE}}, ...}
  \item{cdf}{An \code{\link{AffymetrixCdfFile}} object.}
  \item{...}{Not used.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, already imported arrays will be skipped.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link{CnChipEffectSet}}.
}

\details{
  This import method is robust and memory efficient.  One array at the
  time is imported by first writing to a temporary file which is then
  renamed to the final name, if import was successful.  (If the import
  failed, a temporary file will rename that has to be deleted manually).

  Since only one array at the time is imported, the memory overhead
  will be bounded allowing to import very large tab-delimited data files
  containing a large number of arrays.  Unfortunately, this approach
  slows down the reading substantially, because in each import all but
  one column is parsed but ignored.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{CnChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
