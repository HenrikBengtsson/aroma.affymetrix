%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenericReporter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{GenericReporter}
\docType{class}
\alias{GenericReporter}


\title{The GenericReporter class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class GenericReporter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{GenericReporter}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffymetrixCelSetReporter}}, \emph{\link[aroma.affymetrix]{AffymetrixFileSetReporter}}, \emph{\link[aroma.affymetrix]{SpatialReporter}}\cr

public abstract static class \bold{GenericReporter}\cr
extends \link[R.oo]{Object}\cr


}

\usage{
GenericReporter(tags="*", ...)
}

\arguments{
  \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be added to the output path.}
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getFullName} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getName.GenericReporter]{getName}} \tab Gets the name of the explorer.\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getTags.GenericReporter]{getTags}} \tab Gets the tags of the reporter.\cr
 \tab \code{\link[aroma.affymetrix:process.GenericReporter]{process}} \tab Generates report.\cr
 \tab \code{setup} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
