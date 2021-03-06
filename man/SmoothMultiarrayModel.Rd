%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SmoothMultiarrayModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SmoothMultiarrayModel}
\docType{class}
\alias{SmoothMultiarrayModel}


\title{The SmoothMultiarrayModel class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class SmoothMultiarrayModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ChromosomalModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{SmoothMultiarrayModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{SmoothRmaModel}}, \emph{\link[aroma.affymetrix]{SmoothSaModel}}\cr

public abstract static class \bold{SmoothMultiarrayModel}\cr
extends \emph{\link[aroma.core]{ChromosomalModel}}\cr



 This abstract class represents a chromosomal smoothing method done
 chromosome by chromosome.
}

\usage{
SmoothMultiarrayModel(..., typoOfWeights=c("none", "1/s2"), bandwidth=10000, tags="*")
}

\arguments{
  \item{...}{Arguments passed to the constructor of \code{\link{ChromosomalModel}}.}
  \item{typoOfWeights}{A \code{\link[base]{character}} string.}
  \item{bandwidth}{A single \code{\link[base]{numeric}} specifying the smoothing bandwidth
    in units of nucleotides.}
  \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be added.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getBandwidth} \tab  -\cr
 \tab \code{getOutputTuple} \tab  -\cr
 \tab \code{setBandwidth} \tab  -\cr
}


 \bold{Methods inherited from ChromosomalModel}:\cr
as.character, fit, getAlias, getAromaGenomeTextFile, getAsteriskTags, getChipType, getChipTypes, getChromosomes, getFullName, getFullNames, getGenome, getGenomeData, getGenomeFile, getListOfAromaUgpFiles, getName, getNames, getParentPath, getPath, getReportPath, getRootPath, getSetTuple, getSets, getTags, indexOf, nbrOfArrays, nbrOfChipTypes, setChromosomes, setGenome, getListOfGenomeInformations, getPcuTheta, getPositionChipTypeUnit

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}

\seealso{
 \code{\link{CopyNumberSegmentationModel}}.
}
\keyword{classes}
