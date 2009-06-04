%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChromosomalModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ChromosomalModel}
\docType{class}
\alias{ChromosomalModel}


 \title{The ChromosomalModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ChromosomalModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{ChromosomalModel}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CbsModel}, \emph{\link[aroma.affymetrix]{CopyNumberChromosomalModel}}, \emph{\link[aroma.affymetrix]{CopyNumberSegmentationModel}}, \link[aroma.affymetrix]{GladModel}, \link[aroma.affymetrix]{HaarSegModel}, \emph{\link[aroma.affymetrix]{RawCopyNumberModel}}, \emph{\link[aroma.affymetrix]{SmoothMultiarrayModel}}, \link[aroma.affymetrix]{SmoothRmaModel}, \link[aroma.affymetrix]{SmoothSaModel}\cr

public abstract static class \bold{ChromosomalModel}\cr
extends \link[R.oo]{Object}\cr



  This \emph{abstract} class represents a chromosomal model.
 }

 \usage{ChromosomalModel(cesTuple=NULL, tags="*", genome="Human", ...)}

 \arguments{
   \item{cesTuple}{A \code{\link{ChipEffectSetTuple}}.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
   \item{genome}{A \code{\link[base]{character}} string specifying what genome is process.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{clearCache} \tab  -\cr
 \tab \code{fit} \tab  -\cr
 \tab \code{getAlias} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getArrays.ChromosomalModel]{getArrays}} \tab Gets the names of the arrays.\cr
 \tab \code{getChipEffectFiles} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChipType.ChromosomalModel]{getChipType}} \tab Gets a label for all chip types merged.\cr
 \tab \code{getChipTypes} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChromosomes.ChromosomalModel]{getChromosomes}} \tab Gets the chromosomes available.\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getFullNames} \tab  -\cr
 \tab \code{getGenome} \tab  -\cr
 \tab \code{getListOfAromaUgpFiles} \tab  -\cr
 \tab \code{getListOfChipEffectSets} \tab  -\cr
 \tab \code{getListOfGenomeInformations} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getNames} \tab  -\cr
 \tab \code{getParentPath} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getTableOfArrays.ChromosomalModel]{getTableOfArrays}} \tab Gets a table of arrays.\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.ChromosomalModel]{nbrOfArrays}} \tab Gets the number of arrays.\cr
 \tab \code{\link[aroma.affymetrix:nbrOfChipTypes.ChromosomalModel]{nbrOfChipTypes}} \tab Gets the number of chip types.\cr
 \tab \code{setAlias} \tab  -\cr
 \tab \code{setGenome} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Requirements}{
   This class requires genome information annotation files for
   every chip type.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
