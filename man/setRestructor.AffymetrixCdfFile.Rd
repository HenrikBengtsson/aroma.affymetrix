%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setRestructor.AffymetrixCdfFile}
\alias{setRestructor.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.setRestructor}
\alias{setRestructor,AffymetrixCdfFile-method}

\title{Specifies a function through which}

\description{
 Specifies a function through which of all or a subset of the units.
}

\usage{
\method{setRestructor}{AffymetrixCdfFile}(this, fcn=NULL, ...)
}

\arguments{
  \item{...}{Additional arguments passed to
     \code{\link[affxparser]{readCdfCellIndices}}.}
}

\value{
 Returns the \code{\link[base]{list}} structure returned by
 \code{\link[affxparser]{readCdfCellIndices}}.
}

\section{Requirements}{
  The reconstructor function \emph{have to}:
 \enumerate{
  \item Accept a \code{\link[base]{list}} structure as its first argument.
  \item Require no other arguments.
  \item Return a \code{\link[base]{list}} structure of identical length as the input one.
        In other words, it must not change the number of units nor
        the order of units.
 }

 The reconstructor function \emph{may}:
 \enumerate{
  \item Rearrange the groups or change the number of groups, for
        instance by utilizing \code{\link[affxparser]{applyCdfGroups}}.
  \item Exclude some cells, for instance by excluding MMs.
 }
}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
