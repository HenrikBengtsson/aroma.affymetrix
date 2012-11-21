%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixFileSet}
\docType{class}
\alias{AffymetrixFileSet}


 \title{The AffymetrixFileSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixFileSet}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AffymetrixFileSet}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixCelSet}, \emph{\link[aroma.affymetrix]{AffymetrixCnChpSet}}, \link[aroma.affymetrix]{ChipEffectSet}, \link[aroma.affymetrix]{CnChipEffectSet}, \emph{\link[aroma.affymetrix]{DChipDcpSet}}, \link[aroma.affymetrix]{ExonChipEffectSet}, \link[aroma.affymetrix]{FirmaSet}, \link[aroma.affymetrix]{ParameterCelSet}, \link[aroma.affymetrix]{QualityAssessmentSet}, \link[aroma.affymetrix]{ResidualSet}, \link[aroma.affymetrix]{SnpChipEffectSet}, \link[aroma.affymetrix]{WeightsSet}\cr

public abstract static class \bold{AffymetrixFileSet}\cr
extends \emph{\link[aroma.core]{AromaPlatformInterface}}\cr



  An AffymetrixFileSet object represents a set of \code{\link{AffymetrixFile}}s
  with \emph{identical} chip types.
 }

 \usage{AffymetrixFileSet(files=NULL, ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{AffymetrixFile}}:s.}
   \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFileSet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:as.AffymetrixFileSet.AffymetrixFileSet]{as.AffymetrixFileSet}} \tab Coerce an object to an AffymetrixFileSet object.\cr
 \tab \code{\link[aroma.affymetrix:byPath.AffymetrixFileSet]{byPath}} \tab Defines an AffymetrixFileSet object by searching for Affymetrix files.\cr
}


 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getAromaUgpFile, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile, isCompatibleWith

 \bold{Methods inherited from AromaMicroarrayDataSet}:\cr
as.AromaMicroarrayDataSetList, as.AromaMicroarrayDataSetTuple, getAromaFullNameTranslatorSet, getAverageFile, getChipType, getDefaultFullName, getPlatform, nbrOfArrays, setAttributesBy, setAttributesBySampleAnnotationFile, setAttributesBySampleAnnotationSet, validate

 \bold{Methods inherited from GenericDataFileSet}:\cr
getFullNameTranslatorSet, getParentName, append, appendFiles, appendFullNamesTranslator, appendFullNamesTranslatorBydata.frame, appendFullNamesTranslatorByfunction, appendFullNamesTranslatorBylist, appendFullNamesTranslatorByNULL, appendFullNamesTranslatorByTabularTextFile, appendFullNamesTranslatorByTabularTextFileSet, as.list, byName, byPath, clearCache, clearFullNamesTranslator, copyTo, equals, extract, findByName, getChecksum, getChecksumObjects, getDefaultFullName, getFile, getFileClass, getFileSize, getFullNames, getNames, getPath, getPathnames, getSubdirs, hasFile, indexOf, lapply, nbrOfFiles, sapply, seq, setFullNamesTranslator, sortBy, update2, updateFullName, updateFullNames, validate

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, print


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
