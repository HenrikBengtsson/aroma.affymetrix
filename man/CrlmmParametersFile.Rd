%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CrlmmParametersFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CrlmmParametersFile}
\docType{class}
\alias{CrlmmParametersFile}


 \title{The CrlmmParametersFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CrlmmParametersFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[R.filesets]{GenericTabularFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaTabularBinaryFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaUnitSignalBinaryFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{CrlmmParametersFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CrlmmParametersFile}\cr
extends \link[aroma.core]{AromaUnitSignalBinaryFile}\cr



  An CrlmmParametersFile is a \code{\link[aroma.core]{AromaUnitSignalBinaryFile}}.
 }

 \usage{CrlmmParametersFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link[aroma.core]{AromaUnitSignalBinaryFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{allocate} \tab  -\cr
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{readParameter} \tab  -\cr
 \tab \code{updateParameter} \tab  -\cr
}


 \bold{Methods inherited from AromaUnitSignalBinaryFile}:\cr
allocate, allocateFromUnitAnnotationDataFile, allocateFromUnitNamesFile, as.character, extractMatrix, extractRawGenomicSignals, fromFile, getAromaUgpFile, getChipType, getExtensionPattern, getFilenameExtension, getPlatform, isAverageFile, nbrOfUnits, readDataFrame, writeDataFrame

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getAromaUgpFile, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile, isCompatibleWith

 \bold{Methods inherited from AromaTabularBinaryFile}:\cr
[, [[, [<-, allocate, as.character, colMeans, colMedians, colStats, colSums, dimnames<-, getBytesPerColumn, getColClasses, getColumnNames, getRootName, importFrom, lapply, nbrOfColumns, nbrOfRows, readColumns, readDataFrame, readFooter, readHeader, readRawFooter, setAttributesByTags, subset, summary, updateData, updateDataColumn, writeFooter, writeRawFooter

 \bold{Methods inherited from GenericTabularFile}:\cr
as.character, dim, extractMatrix, getColumnNames, getColumnNameTranslator, nbrOfColumns, nbrOfRows, readColumns, readDataFrame, setColumnNameTranslator, translateColumnNames

 \bold{Methods inherited from GenericDataFile}:\cr
getParentName, compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFilename, getFilenameExtension, getFileSize, getFileType, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
