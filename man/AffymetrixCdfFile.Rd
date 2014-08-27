%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffymetrixCdfFile}
\docType{class}
\alias{AffymetrixCdfFile}


\title{The AffymetrixCdfFile class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class AffymetrixCdfFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AromaChipTypeAnnotationFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{UnitAnnotationDataFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{UnitTypesFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{UnitNamesFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AffymetrixCdfFile}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{AffymetrixCdfFile}\cr
extends \emph{\link[aroma.core]{UnitNamesFile}}\cr



 An AffymetrixCdfFile object represents a generic Affymetrix CDF file.
}

\usage{
AffymetrixCdfFile(...)
}

\arguments{
  \item{...}{Arguments passed to \code{\link{AromaChipTypeAnnotationFile}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:convert.AffymetrixCdfFile]{convert}} \tab Converts a CDF into the same CDF but with another format.\cr
 \tab \code{\link[aroma.affymetrix:createExonByTranscriptCdf.AffymetrixCdfFile]{createExonByTranscriptCdf}} \tab Creates an exon-by-transcript CDF.\cr
 \tab \code{getACSFile} \tab  -\cr
 \tab \code{getAromaCellSequenceFile} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getCellIndices.AffymetrixCdfFile]{getCellIndices}} \tab Gets the cell indices unit by unit.\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getDimension} \tab  -\cr
 \tab \code{getFileFormat} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getGenomeInformation.AffymetrixCdfFile]{getGenomeInformation}} \tab Gets genome information for this chip type.\cr
 \tab \code{getImage} \tab  -\cr
 \tab \code{getPlatform} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getUnitNames.AffymetrixCdfFile]{getUnitNames}} \tab Gets the names of each unit.\cr
 \tab \code{\link[aroma.affymetrix:getUnitTypes.AffymetrixCdfFile]{getUnitTypes}} \tab Gets the types of a set of units.\cr
 \tab \code{hasUnitTypes} \tab  -\cr
 \tab \code{isMonocellCdf} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:isPm.AffymetrixCdfFile]{isPm}} \tab Checks which cells (probes) are PMs and not.\cr
 \tab \code{isUniqueCdf} \tab  -\cr
 \tab \code{nbrOfCells} \tab  -\cr
 \tab \code{nbrOfColumns} \tab  -\cr
 \tab \code{nbrOfQcUnits} \tab  -\cr
 \tab \code{nbrOfRows} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{readDataFrame} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:readUnits.AffymetrixCdfFile]{readUnits}} \tab Reads CDF data unit by unit.\cr
}


 \bold{Methods inherited from UnitNamesFile}:\cr
getUnitNames, indexOf, nbrOfUnits

 \bold{Methods inherited from UnitTypesFile}:\cr
getUnitTypes, nbrOfUnits

 \bold{Methods inherited from UnitAnnotationDataFile}:\cr
byChipType, getAromaUflFile, getAromaUgpFile, getChipType, getDefaultExtension, getPlatform, nbrOfUnits

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getAromaUflFile, getAromaUgpFile, getChipType, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile, isCompatibleWith

 \bold{Methods inherited from AromaChipTypeAnnotationFile}:\cr
as.character, byChipType, byName, findByChipType, fromFile, getChipType, getDefaultExtension, getHeader, getPlatform

 \bold{Methods inherited from AromaMicroarrayDataFile}:\cr
getAttributeXY, getChipType, getPlatform, getPloidy, getXAM, hasAttributeXY, isAverageFile, setAttributeXY, setAttributesByTags

 \bold{Methods inherited from FileCacheKeyInterface}:\cr
getCacheKey

 \bold{Methods inherited from CacheKeyInterface}:\cr
getCacheKey

 \bold{Methods inherited from GenericDataFile}:\cr
as.character, clone, compareChecksum, copyTo, equals, fromFile, getAttribute, getAttributes, getChecksum, getChecksumFile, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFileSize, getFileType, getFilename, getFilenameExtension, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, is.na, isFile, isGzipped, linkTo, readChecksum, renameTo, renameToUpperCaseExt, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validate, validateChecksum, writeChecksum, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, asThis


}

\author{Henrik Bengtsson, Ken Simpson}
\keyword{classes}
