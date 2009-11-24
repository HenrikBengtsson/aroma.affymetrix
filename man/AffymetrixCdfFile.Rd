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
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AromaChipTypeAnnotationFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{UnitAnnotationDataFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{UnitTypesFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{UnitNamesFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AffymetrixCdfFile}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{AffymetrixCdfFile}\cr
extends \emph{\link[aroma.core]{UnitNamesFile}}\cr



  An AffymetrixCdfFile object represents a generic Affymetrix CDF file.
 }

 \usage{AffymetrixCdfFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AromaChipTypeAnnotationFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:convert.AffymetrixCdfFile]{convert}} \tab Converts a CDF into the same CDF but with another format.\cr
 \tab \code{createMonoCell} \tab  -\cr
 \tab \code{createUnique} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getCellIndices.AffymetrixCdfFile]{getCellIndices}} \tab Gets the cell indices unit by unit.\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getDimension} \tab  -\cr
 \tab \code{getFileFormat} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getGenomeInformation.AffymetrixCdfFile]{getGenomeInformation}} \tab Gets genome information for this chip type.\cr
 \tab \code{getImage} \tab  -\cr
 \tab \code{getMonoCell} \tab  -\cr
 \tab \code{getPlatform} \tab  -\cr
 \tab \code{getUnique} \tab  -\cr
 \tab \code{getUnitGroupNamesFromUgcMap} \tab  -\cr
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
 \tab \code{setGenomeInformation} \tab  -\cr
 \tab \code{setSnpInformation} \tab  -\cr
}


 \bold{Methods inherited from UnitNamesFile}:\cr
getUnitNames, indexOf, nbrOfUnits

 \bold{Methods inherited from UnitTypesFile}:\cr
getUnitTypes, nbrOfUnits

 \bold{Methods inherited from UnitAnnotationDataFile}:\cr
byChipType, getAromaUflFile, getAromaUgpFile, getChipType, getPlatform, nbrOfUnits

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile

 \bold{Methods inherited from AromaChipTypeAnnotationFile}:\cr
byChipType, byName, findByChipType, fromFile, getChipType, getHeader, getPlatform

 \bold{Methods inherited from AromaMicroarrayDataFile}:\cr
getAttributeXY, getChipType, getPlatform, getPloidy, getXAM, hasAttributeXY, isAverageFile, setAttributesByTags, setAttributeXY

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFilename, getFilenameExtension, getFileSize, getFileType, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
