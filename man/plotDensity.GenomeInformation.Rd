%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotDensity.GenomeInformation}
\alias{plotDensity.GenomeInformation}
\alias{GenomeInformation.plotDensity}
\alias{plotDensity.GenomeInformation}
\alias{plotDensity,GenomeInformation-method}


\title{Plots the density of SNPs for a given chromosome}

\description{
  Plots the density of SNPs for a given chromosome.
}

\usage{\method{plotDensity}{GenomeInformation}(this, chromosome, ..., adjust=1/20, xlab=NULL, main=NULL, annotate=TRUE)}

\arguments{
 \item{chromosome}{The chromsome to be displayed.}
 \item{...}{Additional arguments passed to \code{\link[aroma.affymetrix:getPositions.GenomeInformation]{*getPositions}()}.}
 \item{adjust}{A bandwidth parameter for the density estimator.}
 \item{xlab}{The label on the x-axis.  If \code{\link[base]{NULL}}, a default will generated.}
 \item{main}{The title of the plot.  If \code{\link[base]{NULL}}, a default will generated.}
 \item{annotate}{If \code{\link[base:logical]{TRUE}}, the plot is annotated with extra information.}
}

\value{
  Returns (invisibly) the estimated density.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{GenomeInformation}}.
}
\keyword{internal}
\keyword{methods}
