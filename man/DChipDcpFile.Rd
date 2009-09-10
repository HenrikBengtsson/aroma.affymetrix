%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipDcpFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{DChipDcpFile}
\docType{class}
\alias{DChipDcpFile}


 \title{The DChipDcpFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class DChipDcpFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{DChipDcpFile}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{DChipDcpFile}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixFile}}\cr



  A DChipDcpFile object represents a DChip DCP file.
 }

 \usage{DChipDcpFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{dim} \tab  -\cr
 \tab \code{extractTheta} \tab  -\cr
 \tab \code{fromFile} \tab  -\cr
 \tab \code{getCalls} \tab  -\cr
 \tab \code{getExcludes} \tab  -\cr
 \tab \code{getFileFormat} \tab  -\cr
 \tab \code{getHeader} \tab  -\cr
 \tab \code{getNormalizedIntensities} \tab  -\cr
 \tab \code{getRawIntensities} \tab  -\cr
 \tab \code{getThetas} \tab  -\cr
 \tab \code{getThetasAB} \tab  -\cr
 \tab \code{getThetaStds} \tab  -\cr
 \tab \code{hasMbeiData} \tab  -\cr
 \tab \code{hasNormalizedData} \tab  -\cr
 \tab \code{nbrOfCells} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
}


 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getPlatform, getUnitAnnotationDataFile, getUnitNamesFile, getUnitTypesFile

 \bold{Methods inherited from AromaMicroarrayDataFile}:\cr
getAttributeXY, getChipType, getPlatform, getPloidy, hasAttributeXY, setAttributesByTags, setAttributeXY

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getFilename, getFileSize, getFileType, getFullName, getFullNameTranslator, getLastAccessedOn, getLastModifiedOn, getName, getPath, getPathname, getTags, gunzip, gzip, hasBeenModified, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setFullName, setFullNameTranslator, setName, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \seealso{
   \code{\link{DChipDcpSet}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
