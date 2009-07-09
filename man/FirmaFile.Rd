%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FirmaFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{FirmaFile}
\docType{class}
\alias{FirmaFile}


 \title{The FirmaFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class FirmaFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{FirmaFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{FirmaFile}\cr
extends \link[aroma.affymetrix]{ParameterCelFile}\cr



  This class represents scores calculated by the FIRMA algorithm.
 }

 \usage{FirmaFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixCelFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{extractMatrix} \tab  -\cr
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{getCellIndices} \tab  -\cr
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

 \author{Ken Simpson (ksimpson[at]wehi.edu.au).}

 \seealso{
   An object of this class is typically part of a \code{\link{FirmaSet}}.
 }
\keyword{classes}
