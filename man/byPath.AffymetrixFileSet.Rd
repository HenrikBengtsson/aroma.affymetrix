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
 \item{path}{The path where to search for Affymetrix files.}
 \item{pattern}{The filename pattern for match files.
    If \code{\link[base]{NULL}}, filename extensions corresponding to known subclasses
    of the abstract \code{\link{AffymetrixFile}} class are search for.}
 \item{recursive}{If \code{\link[base:logical]{TRUE}}, subdirectories are search recursively,
    otherwise not.}
 \item{...}{Optional arguments passed to the constructor of the
    static (calling) class.}
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
