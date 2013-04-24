%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotMvsA.AffymetrixCelFile}
\alias{plotMvsA.AffymetrixCelFile}
\alias{AffymetrixCelFile.plotMvsA}
\alias{plotMvsA,AffymetrixCelFile-method}

\title{Plots log-ratio versus log-intensity in a scatter plot}

\description{
 Plots log-ratio versus log-intensity in a scatter plot.
}

\usage{
\method{plotMvsA}{AffymetrixCelFile}(this, reference, indices=NULL, pch=176, xlim=c(0, 16), ylim=c(-1, 1) * diff(xlim), xlab=expression(A == 1/2 \%*\% log[2](y[1] * y[2])), ylab=expression(M == log[2](y[1]/y[2])), ..., annotate=TRUE)
}

\arguments{
  \item{reference}{The reference channel, i.e. the denominator in the
    log ratios.}
  \item{indices}{Indices of the probe signals to be plotted.}
  \item{pch}{The plot symbol.}
  \item{xlim,ylim}{The range of the x and the y axes.}
  \item{xlab,ylab}{The labels on the x and the y axes.}
  \item{...}{Additional arguments passed to \code{\link[graphics]{plot}}.}
  \item{annotate}{If \code{\link[base:logical]{TRUE}}, the plot is annotated with information about
    the data plotted, otherwise not.}
}

\value{
 Returns (invisibly) a \code{\link[base]{data.frame}} with the plotted data in the
 first two columns.
}

\author{Henrik Bengtsson}

\seealso{
  \code{\link[aroma.affymetrix:smoothScatterMvsA.AffymetrixCelFile]{*smoothScatterMvsA}()}.
  \code{\link[aroma.affymetrix:plotMvsX.AffymetrixCelFile]{*plotMvsX}()}.
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
