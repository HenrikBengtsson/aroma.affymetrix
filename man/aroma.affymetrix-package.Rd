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
  This package implements classes for files and sets of files for various Affymetrix file formats, e.g. AffymetrixCdfFile, AffymetrixCelFile, and AffymetrixCelSet.  These are designed to be memory efficient but still being fast.  The idea is to keep all data on file and only read data into memory when needed.  Clever caching mechanisms are used to minimize the overhead of data IO.  All of the above is hidden in the package API and for the developer (and the end user), the data is queried as if it lives in memory.  With this design it is only the diskspace that limits the number of arrays that can be analyzed.  To install, do: source("http://www.aroma-project.org/hbLite.R"); hbInstall("aroma.affymetrix").

  This package should be considered to be in an alpha or beta phase.
  You should expect the API to be changing over time.
}

\section{Requirements}{
  This package requires several packages from CRAN and Bioconductor.
  Some of its methods also depend on packages from R-forge and our
  local braju.com repository.  We are moving more
  and more of the packages to CRAN and Bioconductor.
}

\section{Installation and updates}{

  The preferred way to install this package is:
  \preformatted{
   source("http://www.aroma-project.org/hbLite.R")
   hbInstall("aroma.affymetrix")
  }
}

\section{To get started}{
  To get started, see the online user guides and the vignettes
  \url{http://www.aroma-project.org/}.
}

\section{How to cite this package}{
 In order to keep improving and providing support for this project,
 please cite references [1], [2], or any applicable publication listed on
 \url{http://aroma-project.org/publications/}, whenever you publish
 work that have been used any of the Aroma Framework.
}

\author{Henrik Bengtsson, James Bullard, Kasper Hansen, Pierre Neuvial, Elizabeth Purdom, Mark Robinson, Ken Simpson}

\section{License}{
  The releases of this package is licensed under
  LGPL version 2.1 or newer.

  The development code of the packages is under a private licence
  (where applicable) and patches sent to the author fall under the
  latter license, but will be, if incorporated, released under the
  "release" license above.
}

\references{
 [1] H. Bengtsson, K. Simpson, J. Bullard, and K. Hansen, \emph{aroma.affymetrix: A generic framework in R for analyzing small to very large Affymetrix data sets in bounded memory}, Tech Report \#745, Department of Statistics, University of California, Berkeley, February 2008.
 \cr
 [2] H. Bengtsson, R. Irizarry, B. Carvalho, and T. Speed, \emph{Estimation and assessment of raw copy numbers at the single locus level}, Bioinformatics, 2008.
 \cr
 [3] H. Bengtsson, \emph{The R.oo package - Object-Oriented Programming with References Using Standard R Code}, In Kurt Hornik, Friedrich Leisch and Achim Zeileis, editors, Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003), March 20-22, Vienna, Austria. \url{http://www.ci.tuwien.ac.at/Conferences/DSC-2003/Proceedings/}
 \cr

 For a complete list, see \url{http://aroma-project.org/publications/}.
}
\keyword{package}
