%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffymetrixFile}
\docType{class}
\alias{AffymetrixFile}


\title{The abstract AffymetrixFile class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class AffymetrixFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AffymetrixFile}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffymetrixCdfFile}}, \emph{\link[aroma.affymetrix]{AffymetrixCelFile}}, \emph{\link[aroma.affymetrix]{AffymetrixCnChpFile}}, \emph{\link[aroma.affymetrix]{AffymetrixPgfFile}}, \emph{\link[aroma.affymetrix]{AffymetrixTsvFile}}, \emph{\link[aroma.affymetrix]{AromaChipTypeAnnotationFile}}, \emph{\link[aroma.affymetrix]{ChipEffectFile}}, \emph{\link[aroma.affymetrix]{CnChipEffectFile}}, \emph{\link[aroma.affymetrix]{CnProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{CnagCfhFile}}, \emph{\link[aroma.affymetrix]{DChipCdfBinFile}}, \emph{\link[aroma.affymetrix]{DChipDcpFile}}, \emph{\link[aroma.affymetrix]{ExonChipEffectFile}}, \emph{\link[aroma.affymetrix]{ExonProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{FirmaFile}}, \emph{\link[aroma.affymetrix]{ParameterCelFile}}, \emph{\link[aroma.affymetrix]{ProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{QualityAssessmentFile}}, \emph{\link[aroma.affymetrix]{ResidualFile}}, \emph{\link[aroma.affymetrix]{SnpChipEffectFile}}, \emph{\link[aroma.affymetrix]{SnpProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{WeightsFile}}\cr

public abstract static class \bold{AffymetrixFile}\cr
extends \emph{\link[aroma.core]{AromaPlatformInterface}}\cr



 An AffymetrixFile object represents a single Affymetrix file,
 e.g. an Affymetrix CEL file or an Affymetrix CDF file.
 Note that this class is abstract and can not be instanciated, but
 instead you have to use one of the subclasses or the generic
 \code{fromFile()} method.
}

\usage{
AffymetrixFile(...)
}

\arguments{
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFile}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getAromaUgpFile, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile, isCompatibleWith

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

\author{Henrik Bengtsson}

\seealso{
  An object of this class is typically part of an \code{\link{AffymetrixFileSet}}.
}
\keyword{classes}
