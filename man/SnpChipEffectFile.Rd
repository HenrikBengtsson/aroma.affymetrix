%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SnpChipEffectFile}
\docType{class}
\alias{SnpChipEffectFile}


 \title{The SnpChipEffectFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class SnpChipEffectFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{CacheKeyInterface}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{FileCacheKeyInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{SnpChipEffectFile}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{CnChipEffectFile}}\cr

public abstract static class \bold{SnpChipEffectFile}\cr
extends \emph{\link[aroma.affymetrix]{ChipEffectFile}}\cr



  This class represents estimates of chip effects in the probe-level models.
 }

 \usage{SnpChipEffectFile(..., mergeStrands=FALSE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ChipEffectFile}}.}
   \item{mergeStrands}{Specifies if the strands are merged or not for these
      estimates.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{extractTheta} \tab  -\cr
 \tab \code{extractTotalAndFracB} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


 \bold{Methods inherited from ChipEffectFile}:\cr
as.character, extractChromosomalDataFrame, extractMatrix, extractTheta, findUnitsTodo, getAM, getAsFullCelFile, getCellIndices, getCellMapForMainCdf, getExpandedCellMap, getParameters, getUnitGroupCellArrayMap, getUnitGroupCellMatrixMap, getXAM, mergeGroups, readUnits, writeAsFullCelFile

 \bold{Methods inherited from ParameterCelFile}:\cr
extractDataFrame, extractMatrix, readUnits

 \bold{Methods inherited from ParametersInterface}:\cr
getParameters, getParametersAsString, getParameterSets

 \bold{Methods inherited from AffymetrixCelFile}:\cr
allocateFromCdf, as.character, clone, createFrom, extractMatrix, fromFile, getAm, getCdf, getExtensionPattern, getFileFormat, getImage, getUnitNamesFile, getUnitTypesFile, highlight, image270, nbrOfCells, plotDensity, plotImage, plotMvsA, plotMvsX, range, setCdf, smoothScatterMvsA, writeImage

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getAromaUgpFile, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile, isCompatibleWith

 \bold{Methods inherited from AromaMicroarrayDataFile}:\cr
getAttributeXY, getChipType, getPlatform, getPloidy, getXAM, hasAttributeXY, isAverageFile, setAttributesByTags, setAttributeXY

 \bold{Methods inherited from FileCacheKeyInterface}:\cr
getCacheKey

 \bold{Methods inherited from CacheKeyInterface}:\cr
getCacheKey

 \bold{Methods inherited from GenericDataFile}:\cr
getParentName, as.character, clone, compareChecksum, copyTo, equals, fromFile, getAttribute, getAttributes, getChecksum, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFilename, getFilenameExtension, getFileSize, getFileType, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, isFile, isGzipped, readChecksum, renameTo, renameToUpperCaseExt, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validate, validateChecksum, writeChecksum

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically part of a \code{\link{SnpChipEffectSet}}.
 }
\keyword{classes}
