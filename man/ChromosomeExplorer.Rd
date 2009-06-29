%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChromosomeExplorer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ChromosomeExplorer}
\docType{class}
\alias{ChromosomeExplorer}


 \title{The ChromosomeExplorer class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ChromosomeExplorer}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{Explorer}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{ChromosomeExplorer}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ChromosomeExplorer}\cr
extends \emph{\link[aroma.core]{Explorer}}\cr


 }

 \usage{ChromosomeExplorer(model=NULL, zooms=2^(0:7), ..., version=c("3"))}

 \arguments{
   \item{model}{A \code{\link{CopyNumberChromosomalModel}} object.}
   \item{zooms}{An positive \code{\link[base]{integer}} \code{\link[base]{vector}} specifying for which zoom
    levels the graphics should be generated.}
   \item{...}{Not used.}
   \item{version}{The version of the Explorer HTML/Javascript generated/used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{display} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChromosomes.ChromosomeExplorer]{getChromosomes}} \tab Gets the chromosomes available.\cr
 \tab \code{getFullNames} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getModel.ChromosomeExplorer]{getModel}} \tab Gets the model.\cr
 \tab \code{getNames} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getZooms} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:process.ChromosomeExplorer]{process}} \tab Generates image files, scripts and dynamic pages for the explorer.\cr
 \tab \code{\link[aroma.affymetrix:setArrays.ChromosomeExplorer]{setArrays}} \tab Sets the arrays.\cr
 \tab \code{setCytoband} \tab  -\cr
 \tab \code{setZooms} \tab  -\cr
}


 \bold{Methods inherited from Explorer}:\cr
addIncludes, addIndexFile, display, getAlias, getArrays, getArraysOfInput, getAsteriskTags, getFullName, getIncludePath, getMainPath, getName, getNameOfInput, getParallelSafe, getPath, getReportPathPattern, getSampleLayerPrefix, getSubname, getTags, getTagsOfInput, getTemplatePath, nbrOfArrays, process, setAlias, setArrays, setParallelSafe, setReportPathPattern, setSubname, setup, splitByReportPathPattern

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Generating PNG images}{
   In order to get better looking graphs, but also to be able to generate
   bitmap images on systems without direct bitmap support, which is the case
   when running R in batch mode or on Unix without X11 support, images are
   created using the \code{\link[R.utils]{png2}} device (a wrapper for
   \code{bitmap()} immitating \code{png()}).  The \code{png()} is only
   used if \code{png2()}, which requires Ghostscript, does not.
   Note, when images are created using \code{png2()}, the images does
   not appear immediately, although the function call is completed,
   so be patient.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
  \code{\link{CopyNumberChromosomalModel}}.
 }
\keyword{classes}
