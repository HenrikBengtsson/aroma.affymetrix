%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSetTuple.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffymetrixCelSetTuple}
\docType{class}
\alias{AffymetrixCelSetTuple}


\title{The AffymetrixCelSetTuple class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class AffymetrixCelSetTuple}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFileSetList}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataSetTuple}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{AffymetrixCelSetTuple}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{ChipEffectSetTuple}, \link[aroma.affymetrix]{CnChipEffectSetTuple}\cr

public static class \bold{AffymetrixCelSetTuple}\cr
extends \emph{\link[aroma.core]{AromaMicroarrayDataSetTuple}}\cr


}

\usage{
AffymetrixCelSetTuple(..., .setClass="AffymetrixCelSet")
}

\arguments{
  \item{...}{Arguments passed to the constructor of
     \code{\link[aroma.core]{AromaMicroarrayDataSetTuple}}.}
  \item{.setClass}{The name of the class of the input set.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byPath} \tab  -\cr
}


 \bold{Methods inherited from AromaMicroarrayDataSetTuple}:\cr
as, as.AromaMicroarrayDataSetTuple, byPath, getAsteriskTags, getChipTypes, getFullNames, getSets, getTags, indexOf, nbrOfChipTypes

 \bold{Methods inherited from GenericDataFileSetList}:\cr
as, as.GenericDataFileSetList, as.character, as.data.frame, as.list, assertDuplicates, clone, extract, getAsteriskTags, getDefaultFullName, getFileList, getFileListClass, getFullNames, getNames, getSet, getSets, getTags, indexOf, length, nbrOfFiles, nbrOfSets, setTags

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
