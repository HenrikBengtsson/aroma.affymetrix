%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CopyNumberChromosomalModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CopyNumberChromosomalModel}
\docType{class}
\alias{CopyNumberChromosomalModel}


 \title{The CopyNumberChromosomalModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CopyNumberChromosomalModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{ChromosomalModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{CopyNumberChromosomalModel}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CbsModel}, \emph{\link[aroma.affymetrix]{CopyNumberSegmentationModel}}, \link[aroma.affymetrix]{GladModel}, \link[aroma.affymetrix]{HaarSegModel}, \emph{\link[aroma.affymetrix]{RawCopyNumberModel}}\cr

public abstract static class \bold{CopyNumberChromosomalModel}\cr
extends \emph{\link[aroma.affymetrix]{ChromosomalModel}}\cr



  This \emph{abstract} class represents a copy-number model.
 }

 \usage{CopyNumberChromosomalModel(cesTuple=NULL, refTuple=NULL, tags="*", genome="Human", ...)}

 \arguments{
   \item{cesTuple}{A \code{\link{ChipEffectSetTuple}}.}
   \item{refTuple}{An optional \code{\link{ChipEffectSetTuple}} for pairwise
     comparisons.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
   \item{genome}{A \code{\link[base]{character}} string specifying what genome is process.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{clearCache} \tab  -\cr
 \tab \code{extractRawCopyNumbers} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fit.CopyNumberChromosomalModel]{fit}} \tab Fits the model.\cr
 \tab \code{getReferenceSetTuple} \tab  -\cr
 \tab \code{isPaired} \tab  -\cr
 \tab \code{plotAxesLayers} \tab  -\cr
 \tab \code{plotChromosomesLayers} \tab  -\cr
 \tab \code{plotCytobandLayers} \tab  -\cr
 \tab \code{plotFitLayers} \tab  -\cr
 \tab \code{plotGridHorizontalLayers} \tab  -\cr
 \tab \code{plotRawCopyNumbers} \tab  -\cr
 \tab \code{plotSampleLayers} \tab  -\cr
}


 \bold{Methods inherited from ChromosomalModel}:\cr
as.character, clearCache, fit, getAlias, getArrays, getAsteriskTags, getChipEffectFiles, getChipType, getChipTypes, getChromosomes, getFullName, getFullNames, getGenome, getGenomeData, getGenomeFile, getListOfAromaUgpFiles, getListOfChipEffectSets, getListOfGenomeInformations, getName, getNames, getParentPath, getPath, getPcuTheta, getPositionChipTypeUnit, getReportPath, getRootPath, getSetTuple, getTableOfArrays, getTags, nbrOfArrays, nbrOfChipTypes, setAlias, setGenome

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \section{Requirements}{
   This class requires genome information annotation files for
   every chip type.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
