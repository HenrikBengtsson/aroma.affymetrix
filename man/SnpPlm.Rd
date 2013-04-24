%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SnpPlm}
\docType{class}
\alias{SnpPlm}


\title{The SnpPlm interface class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class SnpPlm}\cr

\code{\link[R.oo]{Interface}}\cr
\code{~~|}\cr
\code{~~+--}\code{SnpPlm}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffineCnPlm}}, \emph{\link[aroma.affymetrix]{AffineSnpPlm}}, \emph{\link[aroma.affymetrix]{AvgCnPlm}}, \emph{\link[aroma.affymetrix]{AvgSnpPlm}}, \link[aroma.affymetrix]{CnPlm}, \emph{\link[aroma.affymetrix]{HetLogAddCnPlm}}, \emph{\link[aroma.affymetrix]{HetLogAddSnpPlm}}, \emph{\link[aroma.affymetrix]{MbeiCnPlm}}, \emph{\link[aroma.affymetrix]{MbeiSnpPlm}}, \emph{\link[aroma.affymetrix]{RmaCnPlm}}, \emph{\link[aroma.affymetrix]{RmaSnpPlm}}\cr

public class \bold{SnpPlm}\cr
extends \link[R.oo]{Interface}\cr



 An \code{\link[R.oo]{Interface}} implementing methods special for
 \code{\link{ProbeLevelModel}}s specific to SNP arrays.
}

\usage{
SnpPlm(...)
}

\arguments{
  \item{...}{Not used.}
}

\section{Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{getChipEffectSet} \tab  -\cr
 \tab \code{getMergeStrands} \tab  -\cr
 \tab \code{getProbeAffinityFile} \tab  -\cr
 \tab \code{setMergeStrands} \tab  -\cr
}


 \bold{Methods inherited from Interface}:\cr
extend, print, uses


}

\section{Requirements}{
  Classes inheriting from this \code{\link[R.oo]{Interface}} must provide the
  following fields:
  \itemize{
   \item{mergeStrands}{A \code{\link[base]{logical}} value indicating if strands should be
      merged or not.}
  }
}

\examples{
for (zzz in 0) {

# Setup verbose output
verbose <- Arguments$getVerbose(-2)
timestampOn(verbose)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define an example dataset using this path
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find any SNP dataset
path <- NULL
if (is.null(path))
  break

if (!exists("ds")) {
  ds <- AffymetrixCelSet$fromFiles(path)
}
print(ds)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create a set of various PLMs for this dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("models", mode="list")) {
  mergeStrands <- TRUE
  models <- list(
    rma = RmaSnpPlm(ds, mergeStrands=mergeStrands),
    mbei = MbeiSnpPlm(ds, mergeStrands=mergeStrands)
#   affine = AffineSnpPlm(ds, background=FALSE, mergeStrands=mergeStrands)
  )
}
print(models)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For each model, fit a few units
#
# Note, by fitting the same set of units across models, the internal
# caching mechanisms of aroma.affymetrix makes sure that the data is
# only read into memory once. See log for reading speed.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
units <- 55+1:100

for (model in models) {
  ruler(verbose)
  fit(model, units=units, force=TRUE, verbose=verbose)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For each unit, plot the estimated (thetaB,thetaA) for all models
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Should we plot the on the log scale?
log <- TRUE

# Do only user to press ENTER if more than one unit is plotted
opar <- par(ask=(length(units) > 1))

Alab <- expression(theta[A])
Blab <- expression(theta[B])
if (log) {
  lim <- c(6, 16)
} else {
  lim <- c(0, 2^15)
}

# For each unit...
for (unit in units) {
  # For all models...
  for (kk in seq(along=models)) {
    ces <- getChipEffects(models[[kk]])
    ceUnit <- ces[unit]
    snpName <- names(ceUnit)[1]
    theta <- ceUnit[[1]]
    thetaA <- theta[[1]]$theta
    thetaB <- theta[[2]]$theta
    if (log) {
      thetaA <- log(thetaA, base=2)
      thetaB <- log(thetaB, base=2)
    }

    # Create the plot?
    if (kk == 1) {
      plot(NA, xlim=lim, ylim=lim, xlab=Blab, ylab=Alab, main=snpName)
      abline(a=0, b=1, lty=2)
    }

    # Plot the estimated parameters
    points(thetaB, thetaA, col=kk, pch=19)
  }
} # for (unit ...)

# Reset graphical parameter settings
par(opar)

} # for (zzz in 0)
rm(zzz)
}

\author{Henrik Bengtsson}
\keyword{classes}
