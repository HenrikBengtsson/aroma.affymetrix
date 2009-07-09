%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnagCfhSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CnagCfhSet}
\docType{class}
\alias{CnagCfhSet}


 \title{The CnagCfhSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CnagCfhSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CnagCfhSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CnagCfhSet}\cr
extends \link[R.filesets]{GenericDataFileSet}\cr



  An CnagCfhSet object represents a set of CNAG CFH files
  with \emph{identical} chip types.
 }

 \usage{CnagCfhSet(files=NULL, ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{CnagCfhFile}}:s.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{[} \tab  -\cr
 \tab \code{[[} \tab  -\cr
 \tab \code{append} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:as.CnagCfhSet.CnagCfhSet]{as.CnagCfhSet}} \tab Coerce an object to an CnagCfhSet object.\cr
 \tab \code{byName} \tab  -\cr
 \tab \code{findByName} \tab  -\cr
 \tab \code{fromFiles} \tab  -\cr
 \tab \code{getAverage} \tab  -\cr
 \tab \code{getAverageAsinh} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getAverageFile.CnagCfhSet]{getAverageFile}} \tab Calculates the mean and the standard deviation of the cell signal (intensity, standard deviation etc.) across the CFH set.\cr
 \tab \code{getAverageLog} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getCdf.CnagCfhSet]{getCdf}} \tab Gets the CDF structure for this CFH set.\cr
 \tab \code{getData} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getTimestamps} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:isDuplicated.CnagCfhSet]{isDuplicated}} \tab Identifies duplicated CFH files.\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.CnagCfhSet]{nbrOfArrays}} \tab Gets the number of arrays in the file set.\cr
 \tab \code{range} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:setCdf.CnagCfhSet]{setCdf}} \tab Sets the CDF structure for this CFH set.\cr
}


 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as.list, byName, clearCache, copyTo, equals, extract, findByName, fromFiles, getAlias, getFile, getFileClass, getFileSize, getFullName, getFullNames, getFullNameTranslator, getName, getNames, getPath, getPathnames, getTags, hasFile, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setFullName, setFullNamesTranslator, setFullNameTranslator, setName, setTags, update2, validate

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \seealso{
   \code{\link{CnagCfhFile}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
