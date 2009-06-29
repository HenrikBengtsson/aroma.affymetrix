%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  OpticalBackgroundCorrection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{OpticalBackgroundCorrection}
\docType{class}
\alias{OpticalBackgroundCorrection}


\title{The OpticalBackgroundCorrection class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class OpticalBackgroundCorrection}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{OpticalBackgroundCorrection}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{OpticalBackgroundCorrection}\cr
extends \emph{\link[aroma.affymetrix]{BackgroundCorrection}}\cr



 This class represents "optical" background adjustment.
}

\usage{OpticalBackgroundCorrection(..., minimum=1)}

\arguments{
  \item{...}{Arguments passed to the constructor of
    \code{\link{ProbeLevelTransform}}.}
  \item{minimum}{The minimum signal allowed after adjustment.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.OpticalBackgroundCorrection]{process}} \tab Performs background correction.\cr
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

\author{Ken Simpson (ksimpson[at]wehi.edu.au).}
\keyword{classes}
