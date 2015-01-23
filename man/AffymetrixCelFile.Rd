%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffymetrixCelFile}
\docType{class}
\alias{AffymetrixCelFile}


\title{The AffymetrixCelFile class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class AffymetrixCelFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AffymetrixCelFile}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{ChipEffectFile}}, \emph{\link[aroma.affymetrix]{CnChipEffectFile}}, \emph{\link[aroma.affymetrix]{CnProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{ExonChipEffectFile}}, \emph{\link[aroma.affymetrix]{ExonProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{FirmaFile}}, \emph{\link[aroma.affymetrix]{ParameterCelFile}}, \emph{\link[aroma.affymetrix]{ProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{QualityAssessmentFile}}, \emph{\link[aroma.affymetrix]{ResidualFile}}, \emph{\link[aroma.affymetrix]{SnpChipEffectFile}}, \emph{\link[aroma.affymetrix]{SnpProbeAffinityFile}}, \emph{\link[aroma.affymetrix]{WeightsFile}}\cr

public abstract static class \bold{AffymetrixCelFile}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixFile}}\cr



 An AffymetrixCelFile object represents a single Affymetrix CEL file.
}

\usage{
AffymetrixCelFile(..., cdf=NULL)
}

\arguments{
  \item{...}{Arguments passed to \code{\link[aroma.core]{AromaMicroarrayDataFile}}.}
  \item{cdf}{An optional \code{\link{AffymetrixCdfFile}} making it possible to
    override the default CDF file as specified by the CEL file header.
    The requirement is that its number of cells must match that of
    the CEL file.
    If \code{\link[base]{NULL}}, the CDF structure is inferred from the the chip type
    as specified in the CEL file header.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{extractMatrix} \tab  -\cr
 \tab \code{getAm} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getFileFormat} \tab  -\cr
 \tab \code{getImage} \tab  -\cr
 \tab \code{getUnitNamesFile} \tab  -\cr
 \tab \code{getUnitTypesFile} \tab  -\cr
 \tab \code{image270} \tab  -\cr
 \tab \code{nbrOfCells} \tab  -\cr
 \tab \code{plotDensity} \tab  -\cr
 \tab \code{plotImage} \tab  -\cr
 \tab \code{plotMvsA} \tab  -\cr
 \tab \code{plotMvsX} \tab  -\cr
 \tab \code{setCdf} \tab  -\cr
 \tab \code{smoothScatterMvsA} \tab  -\cr
 \tab \code{writeImage} \tab  -\cr
}


 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getAromaUflFile, getAromaUgpFile, getChipType, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile, isCompatibleWith

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

\section{For developers}{
  If you subclass this class, please make sure to query the
  \code{\link{AffymetrixCdfFile}} object (see \code{\link[aroma.affymetrix:getCdf.AffymetrixCelFile]{*getCdf}()}) whenever
  querying CDF information.  Do not use the CDF file inferred from the
  chip type in CEL header, unless you really want it to be hardwired that
  way, otherwise you will break to possibility to override the CDF
  structure.
}

\author{Henrik Bengtsson}

\seealso{
  An object of this class is typically part of an \code{\link{AffymetrixCelSet}}.
}
\keyword{classes}
