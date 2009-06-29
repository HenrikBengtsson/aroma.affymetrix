%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnChipEffectSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CnChipEffectSet}
\docType{class}
\alias{CnChipEffectSet}


 \title{The CnChipEffectSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CnChipEffectSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpChipEffectSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{CnChipEffectSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CnChipEffectSet}\cr
extends \link[aroma.affymetrix]{SnpChipEffectSet}\cr



  This class represents estimates of chip effects in the probe-level models.
 }

 \usage{CnChipEffectSet(..., combineAlleles="byFirstFile")}

 \arguments{
   \item{...}{Arguments passed to \code{\link{SnpChipEffectSet}}.}
   \item{combineAlleles}{A \code{\link[base]{logical}} indicating if the signals from
      allele A and allele B are combined or not.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{exportTotalAndFracB} \tab  -\cr
 \tab \code{extractTheta} \tab  -\cr
 \tab \code{extractTotalAndFreqB} \tab  -\cr
 \tab \code{getAverageFile} \tab  -\cr
 \tab \code{getBaseline} \tab  -\cr
 \tab \code{getCombineAlleles} \tab  -\cr
 \tab \code{setCombineAlleles} \tab  -\cr
 \tab \code{writeWig} \tab  -\cr
}


 \bold{Methods inherited from SnpChipEffectSet}:\cr
exportTotalAndFracB, extractCNT, extractSnpCnvQSet, extractSnpQSet, extractTheta, extractTotalAndFreqB, fromFiles, getAverageFile, getBaseline, getMergeStrands, setMergeStrands, writeCNT

 \bold{Methods inherited from ChipEffectSet}:\cr
boxplotStats, calculateBaseline, calculateFieldBoxplotStats, calculateNuseBoxplotStats, calculateRleBoxplotStats, extractChromosomalDataFrame, extractMatrix, extractTheta, findByName, findUnitsTodo, fromDataSet, fromFiles, getAM, getAsFullCelSet, getAverageFile, getBaseline, getCellIndices, getParameters, getXAM, plotBoxplot, readUnits, updateUnits

 \bold{Methods inherited from ParameterCelSet}:\cr
extractDataFrame, extractMatrix

 \bold{Methods inherited from AffymetrixCelSet}:\cr
[, [[, append, as, as.AffymetrixCelSet, averageQuantile, byName, convertToUnique, extractAffyBatch, extractMatrix, findByName, fromFiles, fromName, getAverage, getAverageAsinh, getAverageFile, getAverageLog, getCdf, getChipType, getData, getFullName, getIntensities, getPlatform, getTimestamps, getUnitGroupCellMap, getUnitIntensities, getUnitNamesFile, isDuplicated, justSNPRMA, nbrOfArrays, normalizeQuantile, plotDensity, range, readUnits, setCdf, update2, writeSgr

 \bold{Methods inherited from AffymetrixFileSet}:\cr
as, as.AffymetrixFileSet, fromFiles

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getPlatform, getUnitNamesFile

 \bold{Methods inherited from AromaMicroarrayDataSet}:\cr
getChipType, getPlatform, nbrOfArrays, setAttributesBy, setAttributesBySampleAnnotationFile, setAttributesBySampleAnnotationSet, validate

 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as.list, byName, clearCache, copyTo, equals, extract, findByName, fromFiles, getAlias, getFile, getFileClass, getFileSize, getFullName, getFullNames, getFullNameTranslator, getName, getNames, getPath, getPathnames, getTags, hasFile, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setFullName, setFullNamesTranslator, setFullNameTranslator, setName, setTags, update2, validate

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
