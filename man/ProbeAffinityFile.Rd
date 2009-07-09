%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeAffinityFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ProbeAffinityFile}
\docType{class}
\alias{ProbeAffinityFile}


 \title{The ProbeAffinityFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ProbeAffinityFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{ProbeAffinityFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CnProbeAffinityFile}, \link[aroma.affymetrix]{ExonProbeAffinityFile}, \link[aroma.affymetrix]{SnpProbeAffinityFile}\cr

public static class \bold{ProbeAffinityFile}\cr
extends \link[aroma.affymetrix]{ParameterCelFile}\cr



  This class represents estimates of probe affinities in probe-level models.
 }

 \usage{ProbeAffinityFile(..., probeModel=c("pm", "mm", "pm-mm", "min1(pm-mm)", "pm+mm"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ParameterCelFile}}.}
   \item{probeModel}{The specific type of probe model.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


 \bold{Methods inherited from ParameterCelFile}:\cr
extractDataFrame, extractMatrix, getParameters, readUnits

 \bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, allocateFromCdf, clone, createFrom, extractMatrix, fromFile, getAm, getCdf, getFileFormat, getImage, getUnitNamesFile, getUnitTypesFile, highlight, image270, nbrOfCells, plotDensity, plotImage, plotMvsA, plotMvsX, setCdf, smoothScatterMvsA, writeImage

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile

 \bold{Methods inherited from AromaMicroarrayDataFile}:\cr
getAttributeXY, getChipType, getPlatform, getPloidy, hasAttributeXY, setAttributesByTags, setAttributeXY

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getFilename, getFileSize, getFileType, getFullName, getFullNameTranslator, getLastAccessedOn, getLastModifiedOn, getName, getPath, getPathname, getTags, gunzip, gzip, hasBeenModified, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setFullName, setFullNameTranslator, setName, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically obtained through the
   \code{getProbeAffinityFile()} method for the \code{\link{ProbeLevelModel}} class.
 }
\keyword{classes}
