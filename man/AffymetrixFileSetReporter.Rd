%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSetReporter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffymetrixFileSetReporter}
\docType{class}
\alias{AffymetrixFileSetReporter}


\title{The AffymetrixFileSetReporter class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class AffymetrixFileSetReporter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericReporter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{AffymetrixFileSetReporter}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffymetrixCelSetReporter}}, \emph{\link[aroma.affymetrix]{SpatialReporter}}\cr

public abstract static class \bold{AffymetrixFileSetReporter}\cr
extends \emph{\link[aroma.affymetrix]{GenericReporter}}\cr


}

\usage{
AffymetrixFileSetReporter(set=NULL, ..., .setClass="AffymetrixFileSet")
}

\arguments{
  \item{set}{An \code{\link{AffymetrixFileSet}} object.}
  \item{...}{Arguments passed to \code{\link{GenericReporter}}.}
  \item{.setClass}{The name of the class of the input set.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getInputName} \tab  -\cr
 \tab \code{getInputTags} \tab  -\cr
}


 \bold{Methods inherited from GenericReporter}:\cr
as.character, getAlias, getAsteriskTags, getFullName, getInputName, getInputTags, getMainPath, getName, getPath, getReportSet, getRootPath, getTags, process, setAlias, setup

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
