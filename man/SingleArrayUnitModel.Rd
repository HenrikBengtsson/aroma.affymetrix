%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SingleArrayUnitModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SingleArrayUnitModel}
\docType{class}
\alias{SingleArrayUnitModel}


 \title{The SingleArrayUnitModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class SingleArrayUnitModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{SingleArrayUnitModel}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{SingleArrayUnitModel}\cr
extends \emph{\link[aroma.affymetrix]{UnitModel}}\cr



  This abstract class represents a unit model that fits one model per unit
  based on signals from a single arrays.
  The nature of a single-array unit model is that each array can be fitted
  independently of the others.
 }

 \usage{SingleArrayUnitModel(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{UnitModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:fit.SingleArrayUnitModel]{fit}} \tab Estimates the model parameters.\cr
}


 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getAsteriskTags, getFitSingleCellUnitFunction

 \bold{Methods inherited from Model}:\cr
fit, getAlias, getAsteriskTags, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
