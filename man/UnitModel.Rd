%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  UnitModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{UnitModel}
\docType{class}
\alias{UnitModel}


 \title{The UnitModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class UnitModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{UnitModel}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffineCnPlm}, \link[aroma.affymetrix]{AffinePlm}, \link[aroma.affymetrix]{AffineSnpPlm}, \emph{\link[aroma.affymetrix]{AlleleSummation}}, \link[aroma.affymetrix]{AvgCnPlm}, \link[aroma.affymetrix]{AvgPlm}, \link[aroma.affymetrix]{AvgSnpPlm}, \link[aroma.affymetrix]{ExonRmaPlm}, \link[aroma.affymetrix]{FirmaModel}, \link[aroma.affymetrix]{HetLogAddCnPlm}, \link[aroma.affymetrix]{HetLogAddPlm}, \link[aroma.affymetrix]{HetLogAddSnpPlm}, \link[aroma.affymetrix]{MbeiCnPlm}, \link[aroma.affymetrix]{MbeiPlm}, \link[aroma.affymetrix]{MbeiSnpPlm}, \emph{\link[aroma.affymetrix]{MultiArrayUnitModel}}, \link[aroma.affymetrix]{ProbeLevelModel}, \link[aroma.affymetrix]{RmaCnPlm}, \link[aroma.affymetrix]{RmaPlm}, \link[aroma.affymetrix]{RmaSnpPlm}, \emph{\link[aroma.affymetrix]{SingleArrayUnitModel}}\cr

public abstract static class \bold{UnitModel}\cr
extends \emph{\link[aroma.affymetrix]{Model}}\cr



  This class is abstract and represents a generic unit model, i.e.
  a model that is applied to each unit separately.
 }

 \usage{UnitModel(dataSet=NULL, probeModel=c("pm", "mm", "pm-mm", "min1(pm-mm)", "pm+mm"), shift=0, ...)}

 \arguments{
   \item{dataSet}{An \code{\link{AffymetrixCelSet}} to which this model should
      be fitted.}
   \item{probeModel}{A \code{\link[base]{character}} string specifying how PM and MM values
      should be modelled.  By default only PM signals are used.}
   \item{shift}{An optional amount the signals should be shifted
      (translated) before fitting the model.}
   \item{...}{Arguments passed to the constructor of \code{\link{Model}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:findUnitsTodo.UnitModel]{findUnitsTodo}} \tab Identifies non-fitted units.\cr
 \tab \code{getFitSingleCellUnitFunction} \tab  -\cr
}


 \bold{Methods inherited from Model}:\cr
fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
