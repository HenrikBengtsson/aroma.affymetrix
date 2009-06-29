%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaChipTypeAnnotationFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AromaChipTypeAnnotationFile}
\docType{class}
\alias{AromaChipTypeAnnotationFile}


 \title{The AromaChipTypeAnnotationFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AromaChipTypeAnnotationFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AromaChipTypeAnnotationFile}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixCdfFile}\cr

public abstract static class \bold{AromaChipTypeAnnotationFile}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixFile}}\cr



  An AromaChipTypeAnnotationFile object represents an annotation file for a
  specific chip type.
 }

 \usage{AromaChipTypeAnnotationFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:byChipType.AromaChipTypeAnnotationFile]{byChipType}} \tab Defines an AromaChipTypeAnnotationFile object by chip type.\cr
 \tab \code{\link[aroma.affymetrix:findByChipType.AromaChipTypeAnnotationFile]{findByChipType}} \tab Locates an annotation file by its chip type.\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getHeader.AromaChipTypeAnnotationFile]{getHeader}} \tab Gets the header of the annotation file.\cr
 \tab \code{getPlatform} \tab  -\cr
}


 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getPlatform, getUnitNamesFile

 \bold{Methods inherited from AromaMicroarrayDataFile}:\cr
getAttributeXY, getChipType, getPlatform, getPloidy, hasAttributeXY, setAttributesByTags, setAttributeXY

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getFilename, getFileSize, getFileType, getFullName, getFullNameTranslator, getLastAccessedOn, getLastModifiedOn, getName, getPath, getPathname, getTags, gunzip, gzip, hasBeenModified, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setFullName, setFullNameTranslator, setName, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
