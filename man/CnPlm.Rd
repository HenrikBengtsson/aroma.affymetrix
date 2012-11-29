%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CnPlm}
\docType{class}
\alias{CnPlm}


 \title{The CnPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CnPlm}\cr

\code{\link[R.oo]{Interface}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{SnpPlm}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CnPlm}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffineCnPlm}}, \emph{\link[aroma.affymetrix]{AvgCnPlm}}, \emph{\link[aroma.affymetrix]{HetLogAddCnPlm}}, \emph{\link[aroma.affymetrix]{MbeiCnPlm}}, \emph{\link[aroma.affymetrix]{RmaCnPlm}}\cr

public class \bold{CnPlm}\cr
extends \link[aroma.affymetrix]{SnpPlm}\cr



  This suport class represents a \code{\link{SnpPlm}} specially designed for
  copy-number analysis.
 }

 \usage{CnPlm(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{SnpPlm}}.}
 }

 \section{Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{getChipEffectSet} \tab  -\cr
 \tab \code{getCombineAlleles} \tab  -\cr
 \tab \code{getProbeAffinityFile} \tab  -\cr
 \tab \code{setCombineAlleles} \tab  -\cr
}


 \bold{Methods inherited from SnpPlm}:\cr
getCellIndices, getChipEffectSet, getMergeStrands, getParameters, getProbeAffinityFile, setMergeStrands

 \bold{Methods inherited from Interface}:\cr
extend, print, uses


 }

 \details{
   Models implementing this copy-number PLM, provides either
   allele-specific or total copy-number estimates.
   For allele-specific CNs the underlying \code{\link{SnpPlm}} model is fitted as
   is, i.e. for each allele seperately with or without the strands first
   being merged.

   For total CNs the probe signals for the two alleles are combined
   (=summed; not averaged) on the intensity scale before fitting
   underlying \code{\link{SnpPlm}} model, again with or without the strands
   first being merged.
 }

 \section{Requirements}{
   Classes inheriting from this \code{\link[R.oo]{Interface}} must provide the
   following fields, in addition to the ones according to \code{\link{SnpPlm}}:
   \itemize{
    \item{combineAlleles}{A \code{\link[base]{logical}} indicating if total or allele-specific
      copy numbers should be estimated according to the above averaging.}
   }
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
