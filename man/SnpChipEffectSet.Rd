%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SnpChipEffectSet}
\docType{class}
\alias{SnpChipEffectSet}


\title{The SnpChipEffectSet class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class SnpChipEffectSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataSet}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{SnpChipEffectSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CnChipEffectSet}\cr

public static class \bold{SnpChipEffectSet}\cr
extends \link[aroma.affymetrix]{ChipEffectSet}\cr



 This class represents estimates of chip effects in the probe-level models.
}

\usage{
SnpChipEffectSet(..., mergeStrands="byFirstFile")
}

\arguments{
  \item{...}{Arguments passed to \code{\link{ChipEffectSet}}.}
  \item{mergeStrands}{Specifies if the strands are merged or not for these
     estimates.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{extractAlleleSet} \tab  -\cr
 \tab \code{extractSnpCnvQSet} \tab  -\cr
 \tab \code{extractSnpQSet} \tab  -\cr
 \tab \code{extractTheta} \tab  -\cr
 \tab \code{extractTotalAndFreqB} \tab  -\cr
 \tab \code{getAverageFile} \tab  -\cr
 \tab \code{getMergeStrands} \tab  -\cr
 \tab \code{setMergeStrands} \tab  -\cr
}


 \bold{Methods inherited from ChipEffectSet}:\cr
as.character, boxplotStats, byPath, calculateBaseline, calculateFieldBoxplotStats, calculateNuseBoxplotStats, calculateRleBoxplotStats, extractAffyBatch, extractChromosomalDataFrame, extractExpressionSet, extractMatrix, extractTheta, findByName, findUnitsTodo, fromDataSet, getAM, getAsFullCelSet, getAverageFile, getBaseline, getCellIndices, getXAM, plotBoxplot, readUnits, updateUnits

 \bold{Methods inherited from ParameterCelSet}:\cr
extractDataFrame, extractMatrix

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from AffymetrixCelSet}:\cr
append, as, as.AffymetrixCelSet, as.character, averageQuantile, byName, byPath, clone, convertToUnique, doCRMAv1, doCRMAv2, doFIRMA, doGCRMA, doRMA, extractAffyBatch, extractFeatureSet, extractMatrix, extractSnpFeatureSet, findByName, getAverage, getAverageAsinh, getAverageFile, getAverageLog, getCdf, getChipType, getData, getIntensities, getPlatform, getTimestamps, getUnitGroupCellMap, getUnitIntensities, getUnitNamesFile, getUnitTypesFile, isDuplicated, justRMA, justSNPRMA, nbrOfArrays, normalizeQuantile, plotDensity, range, readUnits, setCdf, update2, writeSgr

 \bold{Methods inherited from AffymetrixFileSet}:\cr
as, as.AffymetrixFileSet, byPath, getDefaultFullName

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getAromaUflFile, getAromaUgpFile, getChipType, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile, isCompatibleWith

 \bold{Methods inherited from AromaMicroarrayDataSet}:\cr
as.AromaMicroarrayDataSetList, as.AromaMicroarrayDataSetTuple, getAromaFullNameTranslatorSet, getAverageFile, getChipType, getDefaultFullName, getPlatform, setAttributesBy, setAttributesBySampleAnnotationFile, setAttributesBySampleAnnotationSet, validate

 \bold{Methods inherited from GenericDataFileSet}:\cr
[, [[, anyDuplicated, anyNA, append, appendFiles, appendFullNamesTranslator, appendFullNamesTranslatorByNULL, appendFullNamesTranslatorByTabularTextFile, appendFullNamesTranslatorByTabularTextFileSet, appendFullNamesTranslatorBydata.frame, appendFullNamesTranslatorByfunction, appendFullNamesTranslatorBylist, as.character, as.list, byName, byPath, c, clearCache, clearFullNamesTranslator, clone, copyTo, dsApply, dsApplyInPairs, duplicated, equals, extract, findByName, findDuplicated, getChecksum, getChecksumFileSet, getChecksumObjects, getDefaultFullName, getFile, getFileClass, getFileSize, getFiles, getFullNames, getNames, getOneFile, getPath, getPathnames, getSubdirs, gunzip, gzip, hasFile, indexOf, is.na, names, nbrOfFiles, rep, resetFullNames, setFullNamesTranslator, sortBy, unique, update2, updateFullName, updateFullNames, validate, getFullNameTranslatorSet, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
