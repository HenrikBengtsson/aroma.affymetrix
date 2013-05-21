%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{writeImage.AffymetrixCelFile}
\alias{writeImage.AffymetrixCelFile}
\alias{AffymetrixCelFile.writeImage}
\alias{writeImage,AffymetrixCelFile-method}

\title{Writes a spatial image of the signals in the CEL file}

\description{
 Writes a spatial image of the signals in the CEL file.
}

\usage{
\method{writeImage}{AffymetrixCelFile}(this, filename=NULL, fullname=NULL, tags=c("*", "sqrt", "gray"),
  imgFormat="png", path=NULL, field=c("intensities", "stdvs", "pixels"), ..., skip=TRUE,
  verbose=FALSE)
}

\arguments{
  \item{filename}{A \code{\link[base]{character}} string specifying the filename of
    the output file.}
  \item{fullname}{A \code{\link[base]{character}} string specifying the full name of
    the output file.}
  \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of optional tags added to the
    already existing tags of the CEL file.}
  \item{imgFormat}{A \code{\link[base]{character}} string specifying the filename extension
    which also defines the image file format.}
  \item{path}{The path where the image file is stored.}
  \item{...}{Arguments passed to \code{\link[aroma.affymetrix:getImage.AffymetrixCelFile]{*getImage}()}.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns the pathname to the image file created.
}

\examples{\dontrun{
  yellow.colors <- function(n) { hsv(h=0.15, v=0:(n-1)/(n-1)) }
  df <- getFile(ds, 1)
  writeImage(df, tags="gray", palette=gray.colors(256), xrange=c(0,200))
  writeImage(df, tags="yellow", palette=yellow.colors(256), xrange=c(0,200))
  writeImage(df, tags="heat", palette=heat.colors(256), xrange=c(0,200))
  writeImage(df, tags="terrain", palette=terrain.colors(256), xrange=c(0,200))
  writeImage(df, tags="topo", palette=topo.colors(256), xrange=c(0,200))
  writeImage(df, tags="cm", palette=cm.colors(256), xrange=c(0,200))
  writeImage(df, tags="rainbow", palette=rainbow(256), xrange=c(0,200))
}}

\author{Ken Simpson, Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
