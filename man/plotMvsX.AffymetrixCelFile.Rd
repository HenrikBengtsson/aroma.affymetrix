%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotMvsX.AffymetrixCelFile}
\alias{plotMvsX.AffymetrixCelFile}
\alias{AffymetrixCelFile.plotMvsX}
\alias{plotMvsX,AffymetrixCelFile-method}

\title{Plots log-ratio versus another variable in a scatter plot}

\description{
 Plots log-ratio versus another variable in a scatter plot.
}

\usage{
\method{plotMvsX}{AffymetrixCelFile}(this, reference, x, indices=NULL, pch=176, ylim=c(-1, 1) * 2, ylab=NULL, ...,
  what=c("M", "A"), add=FALSE, annotate=!add)
}

\arguments{
  \item{reference}{The reference channel, i.e. the denominator in the
    log ratios.}
  \item{x}{The other variable.  A \code{\link[base]{double}} \code{\link[base]{vector}}.}
  \item{indices}{Indices of the probe signals to be plotted.}
  \item{pch}{The plot symbol.}
  \item{xlim,ylim}{The range of the x and the y axes.}
  \item{xlab,ylab}{The labels on the x and the y axes.}
  \item{...}{Additional arguments passed to \code{\link[graphics]{plot}}().}
  \item{annotate}{If \code{\link[base:logical]{TRUE}}, the plot is annotated with information about
    the data plotted, otherwise not.}
}

\value{
 Returns (invisibly) a \code{\link[base]{data.frame}} with the plotted data in the
 first two columns, and remaining data in the following columns.
}

\author{Henrik Bengtsson}

\seealso{
  \code{\link[aroma.affymetrix:plotMvsA.AffymetrixCelFile]{*plotMvsA}()}.
  \code{\link[aroma.affymetrix:smoothScatterMvsA.AffymetrixCelFile]{*smoothScatterMvsA}()}.
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
