%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SnpInformation}
\docType{class}
\alias{SnpInformation}


 \title{The SnpInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class SnpInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{SnpInformation}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{DChipSnpInformation}}, \emph{\link[aroma.affymetrix]{UflSnpInformation}}\cr

public abstract static class \bold{SnpInformation}\cr
extends \emph{\link[R.filesets]{GenericDataFile}}\cr


 }

 \usage{SnpInformation(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:byChipType.SnpInformation]{byChipType}} \tab Static method to define a genome information set by chip type.\cr
 \tab \code{\link[aroma.affymetrix:fromCdf.SnpInformation]{fromCdf}} \tab Static method to define a genome information set from a CDF.\cr
 \tab \code{fromDataSet} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChipType.SnpInformation]{getChipType}} \tab Gets the chip type of this genome information set.\cr
 \tab \code{\link[aroma.affymetrix:getData.SnpInformation]{getData}} \tab Gets all or a subset of the genome information data.\cr
 \tab \code{getFields} \tab  -\cr
 \tab \code{getFragmentLengths} \tab  -\cr
 \tab \code{getFragmentStarts} \tab  -\cr
 \tab \code{getFragmentStops} \tab  -\cr
 \tab \code{nbrOfEnzymes} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{readDataFrame} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFile}:\cr
getParentName, compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFilename, getFilenameExtension, getFileSize, getFileType, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
