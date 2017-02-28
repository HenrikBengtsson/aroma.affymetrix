%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixProbeTabFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffymetrixProbeTabFile}
\docType{class}
\alias{AffymetrixProbeTabFile}


\title{The AffymetrixProbeTabFile class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class AffymetrixProbeTabFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[R.filesets]{ColumnNamesInterface}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[R.filesets]{GenericTabularFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[R.filesets]{TabularTextFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{CacheKeyInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{FileCacheKeyInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AffymetrixProbeTabFile}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{AffymetrixProbeTabFile}\cr
extends \emph{\link[aroma.core]{AromaPlatformInterface}}\cr



 An AffymetrixProbeTabFile represents an interface to query the data
 contained in an Affymetrix probe tab file, e.g.
 \code{Mapping250K_Nsp_probe_tab}.
}

\usage{
AffymetrixProbeTabFile(...)
}

\arguments{
  \item{...}{Arguments passed to the constructor of \code{\link{AffymetrixFile}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byChipType} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getColumnNames} \tab  -\cr
 \tab \code{getData} \tab  -\cr
}


 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getAromaUflFile, getAromaUgpFile, getChipType, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile, isCompatibleWith

 \bold{Methods inherited from FileCacheKeyInterface}:\cr
getCacheKey

 \bold{Methods inherited from CacheKeyInterface}:\cr
getCacheKey

 \bold{Methods inherited from TabularTextFile}:\cr
as.character, getCommentChar, getDefaultColumnClassPatterns, getDefaultColumnClasses, getDefaultColumnNames, getHeader, getReadArguments, hasColumnHeader, nbrOfLines, nbrOfRows, readColumnNames, readColumns, readDataFrame, readLines, readRawHeader, setCommentChar

 \bold{Methods inherited from GenericTabularFile}:\cr
[, as.character, dim, extractMatrix, head, nbrOfColumns, nbrOfRows, readColumns, readDataFrame, tail, writeColumnsToFiles

 \bold{Methods inherited from ColumnNamesInterface}:\cr
appendColumnNamesTranslator, appendColumnNamesTranslatorByNULL, appendColumnNamesTranslatorBycharacter, appendColumnNamesTranslatorByfunction, appendColumnNamesTranslatorBylist, clearColumnNamesTranslator, clearListOfColumnNamesTranslators, getColumnNames, getColumnNamesTranslator, getDefaultColumnNames, getListOfColumnNamesTranslators, nbrOfColumns, setColumnNames, setColumnNamesTranslator, setListOfColumnNamesTranslators, updateColumnNames

 \bold{Methods inherited from GenericDataFile}:\cr
as.character, clone, compareChecksum, copyTo, equals, fromFile, getAttribute, getAttributes, getChecksum, getChecksumFile, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFileSize, getFileType, getFilename, getFilenameExtension, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, is.na, isFile, isGzipped, linkTo, readChecksum, renameTo, renameToUpperCaseExt, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validate, validateChecksum, writeChecksum, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, asThis


}

\section{About probe-tab files}{
 Probe-tab files are provided by Affymetrix and contains information
 on the probes.  Note that not necessarily all probes are represented
 in the files, e.g. typically only PM probes are given and MM are
 left to be inferred from the PMs.

 The below is an extract of the \code{Mapping250K_Nsp_probe_tab} file
 obtained from the Affymetrix website.  Note that columns are separated
 by TABs.
 \preformatted{
  SNP_A-1780270	1633	2398	3	TTGTTAAGCAAGTGAGTTATTTTAT	f	PM	C
  SNP_A-1780270	1633	2399	3	TTGTTAAGCAAGTGACTTATTTTAT	f	PM	G
  SNP_A-1780270	1951	1780	-4	GGATAAAATAAAATAACTCACTTGC	r	PM	C
  ...
  SNP_A-4241299	2553	1658	4	AAACACATTTTTGGGTCGTAAGGAA	f	PM	G
 }
}

\author{Henrik Bengtsson}


\keyword{classes}
\keyword{internal}
