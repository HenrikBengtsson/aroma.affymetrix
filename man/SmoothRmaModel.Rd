%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SmoothSaModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SmoothRmaModel}
\docType{class}
\alias{SmoothRmaModel}


\title{The SmoothRmaModel class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class SmoothRmaModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ChromosomalModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{SmoothMultiarrayModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{SmoothRmaModel}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{SmoothRmaModel}\cr
extends \emph{\link[aroma.affymetrix]{SmoothMultiarrayModel}}\cr



 This class represents the Chromosomal Smoothing Robust Multichip Analysis
 method.
}

\usage{
SmoothRmaModel(...)
}

\arguments{
  \item{...}{Arguments passed to the constructor of
             \code{\link{SmoothMultiarrayModel}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from SmoothMultiarrayModel}:\cr
as.character, createOutputTuple, fitOneChromosome, getAsteriskTags, getBandwidth, getFitUnitGroupFunction, getOutputTuple, getRootPath, setBandwidth

 \bold{Methods inherited from ChromosomalModel}:\cr
as.character, fit, getAlias, getAromaGenomeTextFile, getAsteriskTags, getChipType, getChipTypes, getChromosomes, getFullName, getFullNames, getGenome, getGenomeData, getGenomeFile, getListOfAromaUgpFiles, getName, getNames, getParentPath, getPath, getReportPath, getRootPath, getSetTuple, getSets, getTags, indexOf, nbrOfArrays, nbrOfChipTypes, setChromosomes, setGenome, getListOfGenomeInformations, getPcuTheta, getPositionChipTypeUnit

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
