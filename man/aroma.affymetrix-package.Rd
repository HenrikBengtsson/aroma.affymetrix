%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.package.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{aroma.affymetrix-package}
\alias{aroma.affymetrix-package}
\alias{aroma.affymetrix}
\docType{package}

\title{Package aroma.affymetrix}


 \description{
   This package implements classes for files and sets of files for various Affymetrix file formats, e.g. AffymetrixCdfFile, AffymetrixCelFile, and AffymetrixCelSet.  These are designed to be memory efficient but still being fast.  The idea is to keep all data on file and only read data into memory when needed.  Clever caching mechanisms are used to minimize the overhead of data IO.  All of the above is hidden in the package API and for the developer (and the end user), the data is queried as if it lives in memory.  With this design it is only the diskspace that limits what number of arrays can be analyzed.

   This package should be considered to be in an alpha or beta phase.
   You should expect the API to be changing over time.
 }

 \section{Requirements}{
   This package requires several packages from CRAN, Bioconductor,
   R-forge and our local braju.com repository.  We are moving more
   and more of the package to CRAN and Bioconductor.
 }

 \section{Installation and updates}{

   To install this package, see instructions at
   \url{http://www.braju.com/R/aroma.affymetrix/}.
 }

 \section{To get started}{
   To get started, see the online user guides and the vignettes
   at \url{http://www.braju.com/R/aroma.affymetrix/}.
 }

 \section{How to cite this package}{
   Please cite references [1] and [2] when using this package.
 }

 \author{
   Henrik Bengtsson with great contributions from
   Ken Simpson, Elizabeth Purdom, and Mark Robinson.
 }

 \section{License}{
   The releases of this package is licensed under
   LGPL version 2.1 or newer.

   The development code of the packages is under a private licence
   (where applicable) and patches sent to the author fall under the
   latter license, but will be, if incorporated, released under the
   "release" license above.
 }

 \references{
  Some of the reference below can be found at
  \url{http://www.maths.lth.se/bioinformatics/publications/}.\cr

  [1] H. Bengtsson, K. Simpson, J. Bullard, and K. Hansen, \emph{aroma.affymetrix: A generic framework in R for analyzing small to very large Affymetrix data sets in bounded memory}, Tech Report \#745, Department of Statistics, University of California, Berkeley, February 2008.
 \cr
  [2] H. Bengtsson, R. Irizarry, B. Carvalho, and T. Speed, \emph{Estimation and assessment of raw copy numbers at the single locus level}, Bioinformatics, 2008.
 \cr
  [3] H. Bengtsson, \emph{The R.oo package - Object-Oriented Programming with References Using Standard R Code}, In Kurt Hornik, Friedrich Leisch and Achim Zeileis, editors, Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003), March 20-22, Vienna, Austria. \url{http://www.ci.tuwien.ac.at/Conferences/DSC-2003/Proceedings/}
 \cr
 }
\keyword{package}
