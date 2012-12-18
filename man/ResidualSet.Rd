%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ResidualSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ResidualSet}
\docType{class}
\alias{ResidualSet}


 \title{The ResidualSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ResidualSet}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{ResidualSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ResidualSet}\cr
extends \link[aroma.core]{ParametersInterface}\cr



  This class represents probe-level residuals from probe-level models.
 }

 \usage{ResidualSet(..., probeModel=c("pm"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixCelSet}}.}
   \item{probeModel}{The specific type of model, e.g. \code{"pm"}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{getAverageFile} \tab  -\cr
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


 \bold{Methods inherited from ParametersInterface}:\cr
getParameters, getParametersAsString, getParameterSets

 \bold{Methods inherited from AffymetrixCelSet}:\cr
append, as, as.AffymetrixCelSet, as.character, averageQuantile, byName, byPath, clone, convertToUnique, doCRMAv1, doCRMAv2, doFIRMA, doGCRMA, doRMA, extractAffyBatch, extractMatrix, extractSnpFeatureSet, findByName, getAverage, getAverageAsinh, getAverageFile, getAverageLog, getCdf, getChipType, getData, getIntensities, getPlatform, getTimestamps, getUnitGroupCellMap, getUnitIntensities, getUnitNamesFile, getUnitTypesFile, isDuplicated, justSNPRMA, nbrOfArrays, normalizeQuantile, plotDensity, range, readUnits, setCdf, update2, writeSgr

 \bold{Methods inherited from AffymetrixFileSet}:\cr
as, as.AffymetrixFileSet, byPath, getDefaultFullName

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getAromaUgpFile, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile, isCompatibleWith

 \bold{Methods inherited from AromaMicroarrayDataSet}:\cr
as.AromaMicroarrayDataSetList, as.AromaMicroarrayDataSetTuple, getAromaFullNameTranslatorSet, getAverageFile, getChipType, getDefaultFullName, getPlatform, nbrOfArrays, setAttributesBy, setAttributesBySampleAnnotationFile, setAttributesBySampleAnnotationSet, validate

 \bold{Methods inherited from GenericDataFileSet}:\cr
getFullNameTranslatorSet, getParentName, append, appendFiles, appendFullNamesTranslator, appendFullNamesTranslatorBydata.frame, appendFullNamesTranslatorByfunction, appendFullNamesTranslatorBylist, appendFullNamesTranslatorByNULL, appendFullNamesTranslatorByTabularTextFile, appendFullNamesTranslatorByTabularTextFileSet, as.character, as.list, byName, byPath, clearCache, clearFullNamesTranslator, clone, copyTo, equals, extract, findByName, getChecksum, getChecksumObjects, getDefaultFullName, getFile, getFileClass, getFiles, getFileSize, getFullNames, getNames, getOneFile, getPath, getPathnames, getSubdirs, hasFile, indexOf, nbrOfFiles, resetFullNames, seq, setFullNamesTranslator, sortBy, update2, updateFullName, updateFullNames, validate

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically obtained through the
   \code{getResidualSet()} method for the \code{\link{ProbeLevelModel}} class.
 }
\keyword{classes}
