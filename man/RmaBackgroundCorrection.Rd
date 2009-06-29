%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RmaBackgroundCorrection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{RmaBackgroundCorrection}
\docType{class}
\alias{RmaBackgroundCorrection}


\title{The RmaBackgroundCorrection class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class RmaBackgroundCorrection}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{BackgroundCorrection}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{RmaBackgroundCorrection}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{RmaBackgroundCorrection}\cr
extends \emph{\link[aroma.affymetrix]{BackgroundCorrection}}\cr



 This class represents the RMA background adjustment function.

}

\usage{RmaBackgroundCorrection(..., addJitter=FALSE, jitterSd=0.2)}

\arguments{
  \item{...}{Arguments passed to the constructor of
    \code{\link{BackgroundCorrection}}.}
  \item{addJitter}{If \code{\link[base:logical]{TRUE}}, Zero-mean gaussian noise is added to the
    signals before being background corrected.}
  \item{jitterSd}{Standard deviation of the jitter noise added.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.RmaBackgroundCorrection]{process}} \tab Performs background correction.\cr
}


 \bold{Methods inherited from BackgroundCorrection}:\cr
process

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\section{Jitter noise}{
  The fitting algorithm of the RMA background correction model may not
  converge if there too many small and discrete signals.  To overcome
  this problem, a small amount of noise may be added to the signals
  before fitting the model.  This is an ad hoc solution that seems to
  work.
  However, add Gaussian noise may generate non-positive signals.
}

\details{
  Internally \code{\link[affy]{bg.adjust}} is used to background correct the
  probe signals.  The default is to background correct PM signals only.
}

\author{Ken Simpson (ksimpson[at]wehi.edu.au).}
\keyword{classes}
