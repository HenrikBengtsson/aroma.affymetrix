%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipGenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{DChipGenomeInformation}
\docType{class}
\alias{DChipGenomeInformation}


 \title{The DChipGenomeInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class DChipGenomeInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{GenomeInformation}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{DChipGenomeInformation}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{DChipGenomeInformation}\cr
extends \emph{\link[aroma.affymetrix]{GenomeInformation}}\cr



  This class represents dChip genome information files, which typically
  contains information about chromosomal locations of the units.
 }

 \usage{DChipGenomeInformation(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenomeInformation}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:byChipType.DChipGenomeInformation]{byChipType}} \tab Defines a DChipGenomeInformation object by chip type.\cr
 \tab \code{readDataFrame} \tab  -\cr
}


 \bold{Methods inherited from GenomeInformation}:\cr
as.character, byChipType, fromCdf, fromDataSet, getChipType, getChromosomes, getChromosomeStats, getData, getPositions, getUnitIndices, getUnitsOnChromosome, getUnitsOnChromosomes, isCompatibleWithCdf, nbrOfUnits, plotDensity, readDataFrame, verify

 \bold{Methods inherited from GenericDataFile}:\cr
getParentName, compareChecksum, copyTo, equals, fromFile, getAttribute, getAttributes, getChecksum, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFilename, getFilenameExtension, getFileSize, getFileType, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validate, validateChecksum, writeChecksum

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, print


 }

 \details{
   The dChip genome information files for various chip types can be
   downloaded from \url{http://www.dchip.org/}.  Put each file in a
   directory named identically as the corresponding chip type under the
   \emph{annotations/} directory, e.g.
   \emph{annotations/Mapping50K\_Hind240/50k hind genome info AfAm
   june 05 hg17.xls}.
   Note that dChip changes the filename and file format slightly between
   chip types, but currently the \code{\link[aroma.affymetrix:byChipType.DChipGenomeInformation]{*byChipType}()} basically searches
   for files with names consisting of \code{"genome info"} or
   \code{"genome_info"}.  At least for the most common chip types, there
   is no need to rename the files in order for this class to recognize them.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
