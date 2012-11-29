%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AllelicCrosstalkCalibration.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AllelicCrosstalkCalibration}
\docType{class}
\alias{AllelicCrosstalkCalibration}


 \title{The AllelicCrosstalkCalibration class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AllelicCrosstalkCalibration}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{AllelicCrosstalkCalibration}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AllelicCrosstalkCalibration}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr



  This class represents a calibration function that transforms the
  probe-level signals such that the signals from the two alleles are
  orthogonal.
  The method fits and calibrates PM signals only.  MM signals will not
  affect the model fitting and are unaffected.
 }

 \usage{AllelicCrosstalkCalibration(dataSet=NULL, ..., model=c("asis", "auto", "CRMA", "CRMAv2"), rescaleBy=c("auto", "groups", "all", "none"), targetAvg=c(2200, 2200), subsetToAvg="-XY", mergeShifts=TRUE, B=1, flavor=c("sfit", "expectile"), alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98, lambda=2, pairBy=c("CDF", "sequence"))}

 \arguments{
   \item{dataSet}{An \code{\link{AffymetrixCelSet}}.}
   \item{...}{Arguments passed to the constructor of
     \code{\link{ProbeLevelTransform}}.}
   \item{model}{A \code{\link[base]{character}} string for quickly specifying default
     parameter settings.}
   \item{rescaleBy}{A \code{\link[base]{character}} string specifying what sets of cells
     should be rescaled towards a target average, if any.
     Default is to rescale all cells together.
     If \code{"none"}, no rescaling is done.}
   \item{targetAvg}{The signal(s) that either the average of the sum
     (if one target value) or the average of each of the alleles
     (if two target values) should have after calibration.
     Only used if \code{rescaleBy != "none"}.}
   \item{subsetToAvg}{The indices of the cells (taken as the intersect of
     existing indices) used to calculate average in order to rescale to
     the target average. If \code{\link[base]{NULL}}, all probes are considered.}
   \item{mergeShifts}{If \code{\link[base:logical]{TRUE}}, the shift of the probe sequence
     relative to the SNP position is ignored, otherwise not.}
   \item{B}{An \code{\link[base]{integer}} specifying by how many nucleotides the allelic
     groups should be stratified by. If zero, all SNPs are put in one
     group.}
   \item{flavor}{A \code{\link[base]{character}} string specifying what algorithm is used
     to fit the crosstalk calibration.}
   \item{alpha, q, Q, lambda}{Model fitting parameters.}
   \item{pairBy}{A \code{\link[base]{character}} string specifying how allele probe pairs
     are identified.}
 }

 \section{What probe signals are updated?}{
   Calibration for crosstalk between allele signals applies by definition
   only SNP units.
   Furthermore, it is only SNP units with two or four unit groups that
   are calibrated.  For instance, in at least on custom SNP CDFs we
   know of, there is a small number of SNP units that have six groups.
   \emph{Currently these units are not calibrated (at all).}
   It is only PM probes that will be calibrated.
   Note that, non-calibrated signals will be saved in the output files.
 }

 \section{What probe signals are used to fit model?}{
   All PM probe pairs are used to fit the crosstalk model.
   In the second step where signals are rescaled to a target average,
   it is possible to specify the set of cells that should be included
   when estimating the target average.
 }

 \section{Important about rescaling towards target average}{
   Rescaling each allele-pair group (e.g. AC, AG, AT, CG, CT, GC)
   towards a target average (\code{rescaleBy="groups"})
   \emph{must not} be used for multi-enzyme chip types,
   e.g. GenomeWideSNP\_6.
   If still done, due to confounded effects of non-perfect enzyme
   mixtures etc, there will be a significant bias between raw CNs
   for SNPs and CN probes.
   Instead, for such chip types \emph{all probe signales} should be
   rescale together towards the target average (\code{rescaleBy="all"}).
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.AllelicCrosstalkCalibration]{process}} \tab Calibrates the data set.\cr
}


 \bold{Methods inherited from ProbeLevelTransform}:\cr
getRootPath

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
as.character, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameters, getParametersAsString, getParameterSets

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
