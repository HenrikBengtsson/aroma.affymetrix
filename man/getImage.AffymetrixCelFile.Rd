%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getImage.AffymetrixCelFile}
\alias{getImage.AffymetrixCelFile}
\alias{AffymetrixCelFile.getImage}
\alias{getImage,AffymetrixCelFile-method}

\title{Creates an RGB Image object from a CEL file}

\description{
 Creates an RGB Image object from a CEL file.
}

\usage{
\method{getImage}{AffymetrixCelFile}(this, other=NULL, transforms=list(sqrt), xrange=c(0, Inf), yrange=xrange,
  zrange=c(0, sqrt(2^16)), field=c("intensities", "stdvs", "pixels"), zoom=1, ...,
  readRectFcn=NULL, verbose=FALSE)
}

\arguments{
  \item{other}{An optional \code{\link{AffymetrixCelFile}} of the same chip type,
     that is used for calculating the ratio (non-logged).  Note, to get
     the log-ratios, the \code{log}() function has to be specified as
     the first transform in the \code{\link[base]{list}} of \code{transformations}.}
  \item{xrange, yrange}{\code{\link[base]{vector}}s of length two specifying the
     (x0,x1) and (y0,y1) regions to be extracted.  If \code{\link[base]{NULL}}, the
     complete regions is used.}
  \item{field}{One of the CEL file fields, i.e. \code{"intensities"},
     \code{stdvs}, or \code{pixels}.}
  \item{zoom}{A \code{\link[base]{numeric}} scale factor in (0,+Inf) for resizing the
    imaging. If \code{1}, no resizing is done.}
  \item{palette}{An optional \code{\link[base]{vector}} of color code.}
  \item{...}{Additional arguments passed to
     \code{\link[aroma.affymetrix:readRawDataRectangle.AffymetrixCelFile]{*readRawDataRectangle}()} and more function.}
  \item{readRectFcn}{A \code{\link[base]{function}} taking arguments 'xrange' and 'yrange',
    or \code{\link[base]{NULL}} for the default read function.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns an Image object as defined by the EBImage package.
  If \code{palette==NULL}, the color code is \code{Grayscale}, otherwise
  \code{TrueColor}.
}

\author{Ken Simpson, Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
