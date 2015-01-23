%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipSnpInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{DChipSnpInformation}
\docType{class}
\alias{DChipSnpInformation}


\title{The DChipSnpInformation class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class DChipSnpInformation}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpInformation}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{DChipSnpInformation}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{DChipSnpInformation}\cr
extends \emph{\link[aroma.affymetrix]{SnpInformation}}\cr



 This class represents dChip genome information files, which typically
 contains information on nucleotide sequences and fragment lengths
 of the units.
}

\usage{
DChipSnpInformation(...)
}

\arguments{
  \item{...}{Arguments passed to \code{\link{SnpInformation}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byChipType} \tab  -\cr
 \tab \code{readDataFrame} \tab  -\cr
}


 \bold{Methods inherited from SnpInformation}:\cr
as.character, byChipType, fromCdf, fromDataSet, getChipType, getData, getFields, getFragmentLengths, getFragmentStarts, getFragmentStops, isCompatibleWithCdf, nbrOfEnzymes, nbrOfUnits, readDataFrame, verify

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

\details{
  The dChip SNP information files for various chip types can be
  downloaded from \url{http://www.dchip.org/}.  Put each file in a
  directory named identically as the corresponding chip type under the
  \emph{annotations/} directory, e.g.
  \emph{annotations/Mapping50K\_Hind240/50k hind snp info AfAm
  june 05 hg17.xls}.
  Note that dChip changes the filename and file format slightly between
  chip types, but currently the \code{\link[aroma.affymetrix:byChipType.DChipSnpInformation]{*byChipType}()} basically searches
  for files with names consisting of \code{"snp info"} or
  \code{"snp_info"}.  At least for the most common chip types, there
  is no need to rename the files in order for this class to recognize them.
}

\author{Henrik Bengtsson}
\keyword{classes}
