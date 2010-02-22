%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GcRmaBackgroundCorrection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{GcRmaBackgroundCorrection}
\docType{class}
\alias{GcRmaBackgroundCorrection}


\title{The GcRmaBackgroundCorrection class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class GcRmaBackgroundCorrection}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{GcRmaBackgroundCorrection}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{GcRmaBackgroundCorrection}\cr
extends \emph{\link[aroma.affymetrix]{BackgroundCorrection}}\cr



 This class represents the GCRMA background adjustment function.

}

\usage{GcRmaBackgroundCorrection(..., indicesNegativeControl=NULL, affinities=NULL, type="fullmodel", opticalAdjust=TRUE, gsbAdjust=TRUE, gsbParameters=NULL)}

\arguments{
  \item{...}{Arguments passed to the constructor of
     \code{\link{ProbeLevelTransform}}.}
  \item{indicesNegativeControl, affinities, type, opticalAdjust,
     gsbAdjust, gsbParameters}{Additional arguments passed to the
     internal \code{bgAdjustGcrma()} method.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.GcRmaBackgroundCorrection]{process}} \tab Performs background correction.\cr
}


 \bold{Methods inherited from BackgroundCorrection}:\cr
process

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


}

\author{Ken Simpson (ksimpson[at]wehi.edu.au).}
\keyword{classes}
