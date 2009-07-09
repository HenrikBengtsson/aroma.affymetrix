%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CrlmmParametersSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CrlmmParametersSet}
\docType{class}
\alias{CrlmmParametersSet}


 \title{The CrlmmParametersSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CrlmmParametersSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericTabularFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaTabularBinarySet}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaUnitSignalBinarySet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{CrlmmParametersSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CrlmmParametersSet}\cr
extends \link[aroma.core]{AromaUnitSignalBinarySet}\cr



  An CrlmmParametersSet object represents a set of
  \code{\link{CrlmmParametersFile}}s with \emph{identical} chip types.
 }

 \usage{CrlmmParametersSet(...)}

 \arguments{
   \item{...}{Arguments passed to
     \code{\link[aroma.core]{AromaUnitSignalBinarySet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byName} \tab  -\cr
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{fromFiles} \tab  -\cr
}


 \bold{Methods inherited from AromaUnitSignalBinarySet}:\cr
byName, findByName, getAromaUgpFile, getChipType, getPlatform, validate

 \bold{Methods inherited from GenericTabularFileSet}:\cr
calculateAverageColumnAcrossFiles, extractMatrix

 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as.list, byName, clearCache, copyTo, equals, extract, findByName, fromFiles, getAlias, getFile, getFileClass, getFileSize, getFullName, getFullNames, getFullNameTranslator, getName, getNames, getPath, getPathnames, getTags, hasFile, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setFullName, setFullNamesTranslator, setFullNameTranslator, setName, setTags, update2, validate

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
