%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  LimmaBackgroundCorrection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{LimmaBackgroundCorrection}
\docType{class}
\alias{LimmaBackgroundCorrection}


\title{The LimmaBackgroundCorrection class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class LimmaBackgroundCorrection}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{LimmaBackgroundCorrection}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{NormExpBackgroundCorrection}\cr

public static class \bold{LimmaBackgroundCorrection}\cr
extends \emph{\link[aroma.affymetrix]{BackgroundCorrection}}\cr



 This class represents the various "background" correction methods
 implemented in the \pkg{limma} package.
}

\usage{LimmaBackgroundCorrection(..., args=NULL, addJitter=FALSE, jitterSd=0.2)}

\arguments{
  \item{...}{Arguments passed to the constructor of
    \code{\link{BackgroundCorrection}}.}
  \item{args}{A \code{\link[base]{list}} of additional arguments passed to the
    correction algorithm.}
  \item{addJitter}{If \code{\link[base:logical]{TRUE}}, Zero-mean gaussian noise is added to the
    signals before being background corrected.}
  \item{jitterSd}{Standard deviation of the jitter noise added.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.LimmaBackgroundCorrection]{process}} \tab Performs background correction.\cr
}


 \bold{Methods inherited from BackgroundCorrection}:\cr
process

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


}

\section{Jitter noise}{
  The fitting algorithm of the normal+exponentital background correction
  model may not converge if there too many small and discrete signals.
  To overcome this problem, a small amount of noise may be added to the
  signals before fitting the model.  This is an ad hoc solution that
  seems to work.
  However, adding Gaussian noise may generate non-positive signals.
}

\details{
  By default, only PM signals are background corrected and MMs are
  left unchanged.
}

\author{Henrik Bengtsson.
        Adopted from RmaBackgroundCorrection by Ken Simpson.}

\seealso{
  Internally, \code{\link[limma]{backgroundCorrect}} is used.
}
\keyword{classes}
