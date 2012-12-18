%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenomeInformation}
\docType{class}
\alias{GenomeInformation}


 \title{The GenomeInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GenomeInformation}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{GenomeInformation}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffymetrixCsvGenomeInformation}}, \emph{\link[aroma.affymetrix]{DChipGenomeInformation}}, \emph{\link[aroma.affymetrix]{UgpGenomeInformation}}\cr

public abstract static class \bold{GenomeInformation}\cr
extends \link[aroma.core]{FileCacheKeyInterface}\cr


 }

 \usage{GenomeInformation(..., .verify=TRUE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFile}}.}
   \item{.verify}{For internal use only.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:byChipType.GenomeInformation]{byChipType}} \tab Static method to define a genome information set by chip type.\cr
 \tab \code{\link[aroma.affymetrix:getChipType.GenomeInformation]{getChipType}} \tab Gets the chip type of this genome information set.\cr
 \tab \code{getChromosomes} \tab  -\cr
 \tab \code{getChromosomeStats} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getData.GenomeInformation]{getData}} \tab Gets all or a subset of the genome information data.\cr
 \tab \code{\link[aroma.affymetrix:getPositions.GenomeInformation]{getPositions}} \tab Gets the physical positions for a set of units.\cr
 \tab \code{getUnitsOnChromosome} \tab  -\cr
 \tab \code{getUnitsOnChromosomes} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:plotDensity.GenomeInformation]{plotDensity}} \tab Plots the density of SNPs for a given chromosome.\cr
}


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
\keyword{classes}
