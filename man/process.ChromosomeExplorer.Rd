%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChromosomeExplorer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{process.ChromosomeExplorer}
\alias{process.ChromosomeExplorer}
\alias{ChromosomeExplorer.process}
\alias{process.ChromosomeExplorer}
\alias{process,ChromosomeExplorer-method}


\title{Generates image files, scripts and dynamic pages for the explorer}

\description{
 Generates image files, scripts and dynamic pages for the explorer.
}

\usage{\method{process}{ChromosomeExplorer}(this, arrays=NULL, chromosomes=NULL, ..., zooms=getZooms(this), layers=FALSE, verbose=FALSE)}

\arguments{
  \item{arrays}{A \code{\link[base]{vector}} of arrays specifying which arrays to
   be considered.  If \code{\link[base]{NULL}}, all are processed.}
  \item{chromosome}{A \code{\link[base]{vector}} of chromosomes specifying which
    chromosomes to be considered.  If \code{\link[base]{NULL}}, all are processed.}
  \item{...}{Not used.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns nothing.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ChromosomeExplorer}}.
}
\keyword{internal}
\keyword{methods}
