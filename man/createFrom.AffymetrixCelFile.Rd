%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.createFrom.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{createFrom.AffymetrixCelFile}
\alias{createFrom.AffymetrixCelFile}
\alias{AffymetrixCelFile.createFrom}
\alias{createFrom,AffymetrixCelFile-method}

\title{Creates a CEL file using another as a template}

\description{
 Creates a CEL file using another as a template.
}

\usage{
\method{createFrom}{AffymetrixCelFile}(this, filename, path=NULL, overwrite=FALSE, skip=!overwrite,
  version=c("4", "3"), methods=c("copy", "create"), clear=FALSE, defValue=0, ...,
  verbose=FALSE)
}

\arguments{
 \item{filename, path}{The filename and path of to the CEL
    file to be created.}
 \item{version}{The file-format version of the CEL file to be created.}
 \item{methods}{If \code{"copy"}, the new file is created as a copy of the
    template file.  If \code{"create"}, the new file is created from
    scratch from the template file.}
 \item{clear}{If \code{\link[base:logical]{TRUE}}, the fields of the CEL file are cleared (zeroed),
    otherwise they contain the same information as the source file.}
 \item{defValue}{A \code{\link[base]{numeric}} value that cleared/allocated elements have.}
 \item{...}{Not used.}
 \item{verbose}{See "R.utils::Verbose".}
}

\value{
 Returns a \code{\link{AffymetrixCelFile}} reference to the new CEL file.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
