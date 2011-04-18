%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MatSmoothing.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{MatSmoothing}
\docType{class}
\alias{MatSmoothing}


\title{The MatSmoothing class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class MatSmoothing}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{MatSmoothing}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{MatSmoothing}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr



 This class represents a function for smoothing data with a trimmed mean.
}

\usage{MatSmoothing(..., design=NULL, probeWindow=300, nProbes=10, meanTrim=0.1)}

\arguments{
  \item{...}{Arguments passed to \code{\link{ProbeLevelTransform}}.}
  \item{design}{A design \code{\link[base]{matrix}}.}
  \item{probeWindow}{Bandwidth to use.  Effectively the width is
     2*probeWindow since it looks probeWindow bases in either direction.}
  \item{nProbes}{The minimum number of probes to calculate a MAT score for.}
  \item{meanTrim}{The amount of trimming of the mean in [0,0.5].}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.MatSmoothing]{process}} \tab Processes the data set.\cr
}


 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


}

\author{Mark Robinson (mrobinson[at]wehi.edu.au).}
\keyword{classes}
