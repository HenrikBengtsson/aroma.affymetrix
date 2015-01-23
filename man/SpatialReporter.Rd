%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpatialReporter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SpatialReporter}
\docType{class}
\alias{SpatialReporter}


\title{The SpatialReporter class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class SpatialReporter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericReporter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSetReporter}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelSetReporter}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{SpatialReporter}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{SpatialReporter}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixCelSetReporter}}\cr



 A SpatialReporter generates image files of spatial representations of
 cell signals for each of the arrays in the input set.
}

\usage{
SpatialReporter(..., reference=NULL)
}

\arguments{
  \item{...}{Arguments passed to \code{\link{AffymetrixCelSetReporter}}.}
  \item{reference}{An optional reference \code{\link{AffymetrixCelFile}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{addColorMap} \tab  -\cr
 \tab \code{getColorMaps} \tab  -\cr
 \tab \code{plotMargins} \tab  -\cr
 \tab \code{process} \tab  -\cr
 \tab \code{setColorMaps} \tab  -\cr
}


 \bold{Methods inherited from AffymetrixCelSetReporter}:\cr
as.character, getChipType, getDataSet, getPath, nbrOfArrays

 \bold{Methods inherited from AffymetrixFileSetReporter}:\cr
getFileSet, getInputName, getInputTags

 \bold{Methods inherited from GenericReporter}:\cr
as.character, getAlias, getAsteriskTags, getFullName, getInputName, getInputTags, getMainPath, getName, getPath, getReportSet, getRootPath, getTags, process, setAlias, setup

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
