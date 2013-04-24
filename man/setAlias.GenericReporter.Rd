%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenericReporter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setAlias.GenericReporter}
\alias{setAlias.GenericReporter}
\alias{GenericReporter.setAlias}
\alias{setAlias,GenericReporter-method}

\title{Sets the alias of the report}

\description{
  Sets the alias of the report.
  If specified, the alias overrides the report name, which is used by
  default.
}

\usage{
\method{setAlias}{GenericReporter}(this, alias=NULL, ...)
}

\arguments{
 \item{alias}{A \code{\link[base]{character}} string for the new alias of the report.
  The alias must consists of valid filename characters, and must not
  contain commas, which are used to separate tags.}
 \item{...}{Not used.}
}

\value{
  Returns nothing.
}

\seealso{
  \code{\link[aroma.affymetrix:getAlias.GenericReporter]{*getAlias}()}.
  \code{\link[aroma.affymetrix:getName.GenericReporter]{*getName}()}.
  For more information see \code{\link{GenericReporter}}.
}
\keyword{internal}
\keyword{methods}
