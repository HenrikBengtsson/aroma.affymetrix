%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParameterCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ParameterCelFile}
\docType{class}
\alias{ParameterCelFile}


 \title{The ParameterCelFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ParameterCelFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{ParameterCelFile}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{ChipEffectFile}}, \emph{\link[aroma.affymetrix]{CnChipEffectFile}}, \emph{\link[aroma.affymetrix]{CnProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{ExonChipEffectFile}}, \emph{\link[aroma.affymetrix]{ExonProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{FirmaFile}}, \emph{\link[aroma.affymetrix]{ProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{ResidualFile}}, \emph{\link[aroma.affymetrix]{SnpChipEffectFile}}, \emph{\link[aroma.affymetrix]{SnpProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{WeightsFile}}\cr

public abstract static class \bold{ParameterCelFile}\cr
extends \link[aroma.core]{ParametersInterface}\cr



  A ParameterCelFile object represents parameter estimates.
 }

 \usage{ParameterCelFile(..., encodeFunction=NULL, decodeFunction=NULL)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixCelFile}}.}
   \item{encodeFunction}{A \code{\link[base]{function}} taking a single \code{\link[base]{list}} structure
      as its argument.}
   \item{decodeFunction}{A \code{\link[base]{function}} taking a single \code{\link[base]{list}} structure
      as its argument.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{extractDataFrame} \tab  -\cr
 \tab \code{extractMatrix} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


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

 \section{File format}{
   The idea behind this class is store data fields which by nature have
   one value per probe (per field) in CEL files.  A perfect example is to
   store probe-affinity estimates and their standard deviations.  There
   is one probe affinity per probe so the structure of a CEL file (and
   its coupled CDF file) is well suited to read/write such information.
 
   Consider a unit group with L probes.  A CEL file stores
   \code{intensities} (L floats), \code{stdvs} (L floats), and
   \code{pixels} (L integers).  Thus, for each probe l=1,...,L, a
   (float, float, integer) tuple is stored.  We can use this for any
   information we want.  If we want a slightly different structure,
   we can choose to encode/decode our structure/information to fit the
   structure of the CEL file.  This abstract class provides transparent
   methods for encoding and decoding such information through methods
   \code{encodeUnitGroup()} and \code{decodeUnitGroup()}.
   By subclassing you can implement different types of data structures.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{classes}
\keyword{IO}
