%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GladModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GladModel}
\docType{class}
\alias{GladModel}


 \title{The GladModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GladModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{ChromosomalModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{CopyNumberChromosomalModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{CopyNumberSegmentationModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{GladModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{GladModel}\cr
extends \emph{\link[aroma.affymetrix]{CopyNumberSegmentationModel}}\cr



  This class represents the Gain and Loss Analysis of DNA regions
  (GLAD) model [1].
  This class can model chip-effect estimates obtained from multiple
  chip types, and not all samples have to be available on all chip types.
 }

 \usage{GladModel(cesTuple=NULL, ...)}

 \arguments{
   \item{cesTuple}{A \code{\link{ChipEffectSetTuple}}.}
   \item{...}{Arguments passed to the constructor of
              \code{\link{CopyNumberSegmentationModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getFitFunction} \tab  -\cr
 \tab \code{writeRegions} \tab  -\cr
}


 \bold{Methods inherited from CopyNumberSegmentationModel}:\cr
fit, getAsteriskTags, getFitFunction, getPositionChipTypeUnit, getRegions, getTags, migrateTool, plot, plotCopyNumberRegionLayers, writeRegions

 \bold{Methods inherited from CopyNumberChromosomalModel}:\cr
as.character, calculateChromosomeStatistics, clearCache, estimateSds, extractRawCopyNumbers, fit, getMatrixChipEffectFiles, getRawCnData, getReferenceSetTuple, getRefSetTuple, isPaired, newPlot, plotAxesLayers, plotChromosomesLayers, plotCytobandLayers, plotFitLayers, plotGridHorizontalLayers, plotRawCopyNumbers, plotSampleLayers

 \bold{Methods inherited from ChromosomalModel}:\cr
as.character, clearCache, fit, getAlias, getArrays, getAsteriskTags, getChipEffectFiles, getChipType, getChipTypes, getChromosomes, getFullName, getFullNames, getGenome, getGenomeData, getGenomeFile, getListOfAromaUgpFiles, getListOfChipEffectSets, getListOfGenomeInformations, getName, getNames, getParentPath, getPath, getPcuTheta, getPositionChipTypeUnit, getReportPath, getRootPath, getSetTuple, getTableOfArrays, getTags, nbrOfArrays, nbrOfChipTypes, setAlias, setGenome

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \details{
   Data from multiple chip types are combined "as is".  This is based
   on the assumption that the relative chip effect estimates are
   non-biased (or at the equally biased across chip types).
   Note that in GLAD there is no way to down weight certain data points,
   which is why we can control for differences in variance across
   chip types.
 }

 \section{Benchmarking}{
   In high-density copy numbers analysis, the most time consuming step
   is fitting the GLAD model.  The complexity of the model grows
   more than linearly (squared? exponentially?) with the number of data
   points in the chromosome and sample being fitted.  This is why it
   take much more than twice the time to fit two chip types together
   than separately.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
  \code{\link{CopyNumberSegmentationModel}}.
 }

 \references{
  [1] Hupe P et al. \emph{Analysis of array CGH data: from signal ratio to
      gain and loss of DNA regions}. Bioinformatics, 2004, 20, 3413-3422.\cr
 }
\keyword{classes}
