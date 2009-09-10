%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenomeInformation}
\docType{class}
\alias{GenomeInformation}


 \title{The GenomeInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GenomeInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{GenomeInformation}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixCsvGenomeInformation}, \link[aroma.affymetrix]{DChipGenomeInformation}, \link[aroma.affymetrix]{UgpGenomeInformation}\cr

public abstract static class \bold{GenomeInformation}\cr
extends \link[R.filesets]{GenericDataFile}\cr


 }

 \usage{GenomeInformation(..., .verify=TRUE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFile}}.}
   \item{.verify}{For internal use only.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:byChipType.GenomeInformation]{byChipType}} \tab Static method to define a genome information set by chip type.\cr
 \tab \code{\link[aroma.affymetrix:fromCdf.GenomeInformation]{fromCdf}} \tab Static method to define a genome information set from a CDF.\cr
 \tab \code{fromDataSet} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChipType.GenomeInformation]{getChipType}} \tab Gets the chip type of this genome information set.\cr
 \tab \code{getChromosomes} \tab  -\cr
 \tab \code{getChromosomeStats} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getData.GenomeInformation]{getData}} \tab Gets all or a subset of the genome information data.\cr
 \tab \code{\link[aroma.affymetrix:getPositions.GenomeInformation]{getPositions}} \tab Gets the physical positions for a set of units.\cr
 \tab \code{getUnitsOnChromosome} \tab  -\cr
 \tab \code{getUnitsOnChromosomes} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:plotDensity.GenomeInformation]{plotDensity}} \tab Plots the density of SNPs for a given chromosome.\cr
}


 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getFilename, getFileSize, getFileType, getFullName, getFullNameTranslator, getLastAccessedOn, getLastModifiedOn, getName, getPath, getPathname, getTags, gunzip, gzip, hasBeenModified, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setFullName, setFullNameTranslator, setName, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
