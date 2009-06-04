%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CopyNumberSegmentationModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CopyNumberSegmentationModel}
\docType{class}
\alias{CopyNumberSegmentationModel}


 \title{The CopyNumberSegmentationModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CopyNumberSegmentationModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{ChromosomalModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{CopyNumberChromosomalModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{CopyNumberSegmentationModel}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CbsModel}, \link[aroma.affymetrix]{GladModel}, \link[aroma.affymetrix]{HaarSegModel}\cr

public abstract static class \bold{CopyNumberSegmentationModel}\cr
extends \emph{\link[aroma.affymetrix]{CopyNumberChromosomalModel}}\cr



  This \emph{abstract} class represents a copy-number segmentation model.
 }

 \usage{CopyNumberSegmentationModel(...)}

 \arguments{
   \item{...}{Arguments passed to constructor
      \code{\link{CopyNumberChromosomalModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:fit.CopyNumberSegmentationModel]{fit}} \tab Fits the model.\cr
 \tab \code{getRegions} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{migrateTool} \tab  -\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{plotCopyNumberRegionLayers} \tab  -\cr
 \tab \code{writeRegions} \tab  -\cr
}


 \bold{Methods inherited from CopyNumberChromosomalModel}:\cr
as.character, calculateChromosomeStatistics, clearCache, estimateSds, extractRawCopyNumbers, fit, getMatrixChipEffectFiles, getRawCnData, getReferenceSetTuple, getRefSetTuple, isPaired, newPlot, plotAxesLayers, plotChromosomesLayers, plotCytobandLayers, plotFitLayers, plotGridHorizontalLayers, plotRawCopyNumbers, plotSampleLayers

 \bold{Methods inherited from ChromosomalModel}:\cr
as.character, clearCache, fit, getAlias, getArrays, getAsteriskTags, getChipEffectFiles, getChipType, getChipTypes, getChromosomes, getFullName, getFullNames, getGenome, getGenomeData, getGenomeFile, getListOfAromaUgpFiles, getListOfChipEffectSets, getListOfGenomeInformations, getName, getNames, getParentPath, getPath, getPcuTheta, getPositionChipTypeUnit, getReportPath, getRootPath, getSetTuple, getTableOfArrays, getTags, nbrOfArrays, nbrOfChipTypes, setAlias, setGenome

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
