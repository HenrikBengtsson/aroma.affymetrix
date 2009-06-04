%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getTimestamp.AffymetrixCelFile}
\alias{getTimestamp.AffymetrixCelFile}
\alias{AffymetrixCelFile.getTimestamp}
\alias{getTimestamp.AffymetrixCelFile}
\alias{getTimestamp,AffymetrixCelFile-method}


\title{Gets the timestamp in the CEL header}

\description{
 Gets the timestamp in the CEL header.
}

\usage{\method{getTimestamp}{AffymetrixCelFile}(this, format="\%m/\%d/\%y \%H:\%M:\%S", ...)}

\arguments{
  \item{format}{The default format string for parsing the time stamp.}
  \item{...}{Not used.}
}

\value{
 Returns a \code{POSIXct} object.
 The parsed string containing the timestamp is returned as
 attribute \code{text}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, \code{\link[base]{strptime}}() is used to parse the time stamp.
  \code{\link[base]{DateTimeClasses}} for more information on \code{POSIXct}.
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
