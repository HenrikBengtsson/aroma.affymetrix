%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  NormExpBackgroundCorrection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{NormExpBackgroundCorrection}
\docType{class}
\alias{NormExpBackgroundCorrection}


\title{The NormExpBackgroundCorrection class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class NormExpBackgroundCorrection}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{LimmaBackgroundCorrection}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{NormExpBackgroundCorrection}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{NormExpBackgroundCorrection}\cr
extends \link[aroma.affymetrix]{LimmaBackgroundCorrection}\cr



 This class represents the normal exponential background correction model.
 Estimators of the \pkg{limma} package is used.
}

\usage{NormExpBackgroundCorrection(..., method=c("rma", "saddle", "mle"))}

\arguments{
  \item{...}{Arguments passed to the constructor of
    \code{\link{LimmaBackgroundCorrection}}.}
  \item{method}{The estimator used, cf. argument \code{normexp.method}
    of \code{\link[limma]{backgroundCorrect}} in \pkg{limma} for more details.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from LimmaBackgroundCorrection}:\cr
process

 \bold{Methods inherited from BackgroundCorrection}:\cr
process

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputDataSetOLD20090509, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\author{Henrik Bengtsson}

\seealso{
  Internally, \code{\link[limma]{backgroundCorrect}} is used.
}
\keyword{classes}
