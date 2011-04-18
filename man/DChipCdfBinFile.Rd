%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipCdfBinFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{DChipCdfBinFile}
\docType{class}
\alias{DChipCdfBinFile}


 \title{The DChipCdfBinFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class DChipCdfBinFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{UnitAnnotationDataFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{UnitNamesFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{DChipCdfBinFile}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{DChipCdfBinFile}\cr
extends \emph{\link[aroma.core]{UnitNamesFile}}\cr



  A DChipCdfBinFile object represents a DChip CDF.bin file.
 }

 \usage{DChipCdfBinFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byChipType} \tab  -\cr
 \tab \code{findByChipType} \tab  -\cr
 \tab \code{fromFile} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getFileFormat} \tab  -\cr
 \tab \code{getHeader} \tab  -\cr
 \tab \code{getPlatform} \tab  -\cr
 \tab \code{getUnitNames} \tab  -\cr
 \tab \code{getUnitSizes} \tab  -\cr
 \tab \code{nbrOfCells} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{readDataFrame} \tab  -\cr
}


 \bold{Methods inherited from UnitNamesFile}:\cr
getUnitNames, indexOf, nbrOfUnits

 \bold{Methods inherited from UnitAnnotationDataFile}:\cr
byChipType, getAromaUflFile, getAromaUgpFile, getChipType, getPlatform, nbrOfUnits

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getAromaUgpFile, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile, isCompatibleWith

 \bold{Methods inherited from AromaMicroarrayDataFile}:\cr
getAttributeXY, getChipType, getPlatform, getPloidy, getXAM, hasAttributeXY, isAverageFile, setAttributesByTags, setAttributeXY

 \bold{Methods inherited from GenericDataFile}:\cr
getParentName, compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFilename, getFilenameExtension, getFileSize, getFileType, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
