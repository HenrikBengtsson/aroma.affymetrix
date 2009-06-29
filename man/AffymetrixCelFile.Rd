%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixCelFile}
\docType{class}
\alias{AffymetrixCelFile}


 \title{The AffymetrixCelFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixCelFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{AffymetrixCelFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{ChipEffectFile}, \link[aroma.affymetrix]{CnChipEffectFile}, \link[aroma.affymetrix]{CnProbeAffinityFile}, \link[aroma.affymetrix]{ExonChipEffectFile}, \link[aroma.affymetrix]{ExonProbeAffinityFile}, \link[aroma.affymetrix]{FirmaFile}, \link[aroma.affymetrix]{ParameterCelFile}, \link[aroma.affymetrix]{ProbeAffinityFile}, \link[aroma.affymetrix]{QualityAssessmentFile}, \link[aroma.affymetrix]{ResidualFile}, \link[aroma.affymetrix]{SnpChipEffectFile}, \link[aroma.affymetrix]{SnpProbeAffinityFile}, \link[aroma.affymetrix]{WeightsFile}\cr

public static class \bold{AffymetrixCelFile}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixFile}}\cr



  An AffymetrixCelFile object represents a single Affymetrix CEL file.
 }

 \usage{AffymetrixCelFile(..., cdf=NULL)}

 \arguments{
   \item{...}{Arguments passed to \code{\link[aroma.core]{AromaMicroarrayDataFile}}.}
   \item{cdf}{An optional \code{\link{AffymetrixCdfFile}} making it possible to
     override the default CDF file as specified by the CEL file header.
     The requirement is that its number of cells must match that of
     the CEL file.
     If \code{\link[base]{NULL}}, the CDF structure is inferred from the the chip type
     as specified in the CEL file header.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{[} \tab  -\cr
 \tab \code{[[} \tab  -\cr
 \tab \code{clone} \tab  -\cr
 \tab \code{extractMatrix} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fromFile.AffymetrixCelFile]{fromFile}} \tab Defines an AffymetrixCelFile object from a CEL file.\cr
 \tab \code{getAm} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getCdf.AffymetrixCelFile]{getCdf}} \tab Gets the CDF structure for this CEL file.\cr
 \tab \code{getFileFormat} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getImage.AffymetrixCelFile]{getImage}} \tab Creates an RGB Image object from a CEL file.\cr
 \tab \code{getUnitNamesFile} \tab  -\cr
 \tab \code{highlight} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:image270.AffymetrixCelFile]{image270}} \tab Displays all or a subset of the data spatially.\cr
 \tab \code{nbrOfCells} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:plotDensity.AffymetrixCelFile]{plotDensity}} \tab Plots the density of the probe signals on the array.\cr
 \tab \code{\link[aroma.affymetrix:plotImage.AffymetrixCelFile]{plotImage}} \tab Displays a spatial plot of a CEL file.\cr
 \tab \code{\link[aroma.affymetrix:plotMvsA.AffymetrixCelFile]{plotMvsA}} \tab Plots log-ratio versus log-intensity in a scatter plot.\cr
 \tab \code{\link[aroma.affymetrix:plotMvsX.AffymetrixCelFile]{plotMvsX}} \tab Plots log-ratio versus another variable in a scatter plot.\cr
 \tab \code{\link[aroma.affymetrix:setCdf.AffymetrixCelFile]{setCdf}} \tab Sets the CDF structure for this CEL file.\cr
 \tab \code{\link[aroma.affymetrix:smoothScatterMvsA.AffymetrixCelFile]{smoothScatterMvsA}} \tab Plots log-ratio versus log-intensity in a smooth scatter plot.\cr
 \tab \code{\link[aroma.affymetrix:writeImage.AffymetrixCelFile]{writeImage}} \tab Writes a spatial image of the signals in the CEL file.\cr
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

 \section{For developers}{
   If you subclass this class, please make sure to query the
   \code{\link{AffymetrixCdfFile}} object (see \code{\link[aroma.affymetrix:getCdf.AffymetrixCelFile]{*getCdf}()}) whenever
   querying CDF information.  Do not use the CDF file inferred from the
   chip type in CEL header, unless you really want it to be hardwired that
   way, otherwise you will break to possibility to override the CDF
   structure.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically part of an \code{\link{AffymetrixCelSet}}.
 }
\keyword{classes}
