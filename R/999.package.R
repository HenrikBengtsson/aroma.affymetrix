#########################################################################/**
# @RdocPackage aroma.affymetrix
#
# \description{
#   @eval "getDescription(aroma.affymetrix)" 
#
#   This package should be considered to be in an alpha or beta phase.
#   You should expect the API to be changing over time.
# }
#
# \section{Requirements}{
#   This package requires quite a few package.  Some like
#   \pkg{R.oo} [1] and \pkg{R.utils} are on CRAN, others like
#   \pkg{affxparser} and \pkg{aroma.light} are on Bioconductor,
#   whereas some like \pkg{R.huge} are only on 
#   \url{http://www.braju.com/R/}.  The reason that this package
#   is not at Bioconductor is that it depends on packages that
#   are in a developmental stage and are therefore not accepted
#   as is on CRAN or Bioconductor.
# }
#
# \section{Installation and updates}{
#
#   To install this package, see instructions at 
#   \url{http://www.braju.com/R/}.
# } 
#
# \section{To get started}{
#   To get started, see:
#   \enumerate{
#     \item @see "SnpPlm" - 
#       Contains an example how to fit different allele-specific 
#       probe-level models (PLMs) to the same data set, and how
#       to display the estimated \eqn{(\theta_B,\theta_A)} "chip effects".
#     \item @see "CnPlm" - 
#       Contains an example how to fit copy-number estimates for
#       a specific chromosome, and how to display the estimates along
#       the chromosome.
#     \item \code{normalizeQuantile()} in @see "AffymetrixCelSet" - 
#       Contains an example how to quantile normalize a set of arrays.
#   }
# }
# 
# \section{Performance}{
#   There is a performance price which we have to pay for not keeping 
#   data in memory but on file.  However, the performance is still
#   quite good, because the underlying read methods provided by
#   the \pkg{affxparser} package that is constantly being optimized
#   for I/O speed.
#
#   Note that it is much faster to access files from a local drive than
#   over a local network.  Thus, you might want to consider to copy
#   files on the network to a temporary local directory and work from
#   there.
# }
#
# \section{How to cite this package}{
#   Please cite references [1] and [2] when using this package.
# }
#
# \section{Wishlist}{
#  Here is a list of features that would be useful, but which I have
#  too little time to add myself. Contributions are appreciated.
#  \itemize{
#    \item At the moment, nothing.
#  }
#
#  If you consider to contribute, make sure it is not already 
#  implemented by downloading the latest "devel" version!
# } 
#
# \author{
#   Henrik Bengtsson with great contributions from
#   Ken Simpson, Elizabeth Purdom and Mark Robinson.
# }
#
# \section{License}{
#   The releases of this package is licensed under 
#   LGPL version 2.1 or newer.
#
#   The development code of the packages is under a private licence 
#   (where applicable) and patches sent to the author fall under the
#   latter license, but will be, if incorporated, released under the
#   "release" license above.
# }
# 
# \references{
#  Some of the reference below can be found at 
#  \url{http://www.maths.lth.se/bioinformatics/publications/}.\cr
#
#  [1] @include "../incl/BengtssonH_etal_2008b.bib.Rdoc" \cr
#  [2] @include "../incl/BengtssonH_etal_2008.bib.Rdoc" \cr
#  [3] @include "../incl/BengtssonH_2003.bib.Rdoc" \cr
# }
#*/#########################################################################

