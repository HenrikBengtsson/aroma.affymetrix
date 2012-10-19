%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{byPath.AffymetrixFileSet}
\alias{byPath.AffymetrixFileSet}
\alias{AffymetrixFileSet.byPath}
\alias{byPath.AffymetrixFileSet}
\alias{byPath,AffymetrixFileSet-method}


\title{Defines an AffymetrixFileSet object by searching for Affymetrix files}

\description{
  Defines an AffymetrixFileSet object by searching for Affymetrix files.
}

\synopsis{byPath.AffymetrixFileSet(static, ..., fileClass="AffymetrixFile")}
\usage{AffymetrixFileSet$byPath(..., fileClass="AffymetrixFile")}

\arguments{
 \item{...}{Additional arguments passed to the constructor of the static
    (calling) class.}
 \item{fileClass}{The name of the \code{\link{GenericDataFile}} class.}
}

\value{
  Returns an \code{\link{AffymetrixFileSet}} object.
}

\section{Reserved filenames}{
  Note that files with names starting with a period \code{.} are not
  searched for.  The reason for this is that such files are reserved for
  internal use of this package.  For instance, the package store average
  signals across CEL files in a file named as \code{.average<something>.CEL}
  in the same directory as the CEL files, and when such a directory is
  scanned we do not want such files to be interpreted as data.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
