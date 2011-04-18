%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  UgpGenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{UgpGenomeInformation}
\docType{class}
\alias{UgpGenomeInformation}


 \title{The UgpGenomeInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class UgpGenomeInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{GenomeInformation}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{UgpGenomeInformation}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{UgpGenomeInformation}\cr
extends \emph{\link[aroma.affymetrix]{GenomeInformation}}\cr



  This class represents Aroma UGP genome information files.
 }

 \usage{UgpGenomeInformation(..., .ugp=NULL, .verify=TRUE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenomeInformation}}.}
   \item{.ugp}{For internal use only.}
   \item{.verify}{For internal use only.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:byChipType.UgpGenomeInformation]{byChipType}} \tab Defines a UgpGenomeInformation object by chip type.\cr
 \tab \code{clearCache} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getChromosomes} \tab  -\cr
 \tab \code{getData} \tab  -\cr
 \tab \code{getUnitsOnChromosome} \tab  -\cr
 \tab \code{isCompatibleWithCdf} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{readDataFrame} \tab  -\cr
}


 \bold{Methods inherited from GenomeInformation}:\cr
byChipType, fromCdf, fromDataSet, getChipType, getChromosomes, getChromosomeStats, getData, getPositions, getUnitIndices, getUnitsOnChromosome, getUnitsOnChromosomes, isCompatibleWithCdf, nbrOfUnits, plotDensity, readDataFrame

 \bold{Methods inherited from GenericDataFile}:\cr
getParentName, compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFilename, getFilenameExtension, getFileSize, getFileType, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
