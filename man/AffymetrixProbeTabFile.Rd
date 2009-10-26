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
\code{~~~~~~~~~~~~+--}\code{\link[R.filesets]{GenericTabularFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[R.filesets]{TabularTextFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AffymetrixProbeTabFile}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{AffymetrixProbeTabFile}\cr
extends \emph{\link[aroma.core]{AromaPlatformInterface}}\cr



  An AffymetrixProbeTabFile represents an interface to query the data
  contained in an Affymetrix probe tab file, e.g.
  \code{Mapping250K_Nsp_probe_tab}.
 }

 \usage{AffymetrixProbeTabFile(...)}

 \arguments{
   \item{...}{Arguments passed to the constructor of \code{\link{AffymetrixFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byChipType} \tab  -\cr
 \tab \code{fromCdf} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getColumnNames} \tab  -\cr
 \tab \code{getData} \tab  -\cr
 \tab \code{translateFullName} \tab  -\cr
}


 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile

 \bold{Methods inherited from TabularTextFile}:\cr
[, as.character, getColumnNames, getHeader, getReadArguments, hasColumnHeader, nbrOfLines, nbrOfRows, readColumnNames, readColumns, readDataFrame, readLines, readRawHeader, writeColumnsToFiles

 \bold{Methods inherited from GenericTabularFile}:\cr
as.character, dim, extractMatrix, getColumnNames, getColumnNameTranslator, nbrOfColumns, nbrOfRows, readColumns, readDataFrame, setColumnNameTranslator, translateColumnNames

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFilename, getFilenameExtension, getFileSize, getFileType, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


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

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{classes}
\keyword{internal}
