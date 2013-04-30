%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.extractFeatureSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{extractFeatureSet.AffymetrixCelSet}
\alias{extractFeatureSet.AffymetrixCelSet}
\alias{AffymetrixCelSet.extractFeatureSet}
\alias{extractFeatureSet,AffymetrixCelSet-method}
\alias{extractFeatureSet}
\alias{AffymetrixCelSet.extractSnpFeatureSet}
\alias{extractSnpFeatureSet.AffymetrixCelSet}
\alias{extractSnpFeatureSet,AffymetrixCelSet-method}

\alias{extractSnpFeatureSet}

\title{Extracts CEL signals an in-memory FeatureSet object}

\description{
 Extracts CEL signals an in-memory FeatureSet object from a \code{\link{AffymetrixCelSet}} object.
 Note that any modifications done to the extracted object will \emph{not}
 be reflected in the original CEL set.
}

\usage{
\method{extractFeatureSet}{AffymetrixCelSet}(this, ..., verbose=FALSE)
}

\arguments{
  \item{...}{Additional argument passed to \code{\link[oligo]{read.celfiles}}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link[oligoClasses]{FeatureSet-class}} object.
}

\author{Henrik Bengtsson}

\seealso{
  Internally \code{\link[oligo]{read.celfiles}} is used to read the data.
  To read \emph{summarized} data as a \code{\link[Biobase]{ExpressionSet-class}}
  object, see \code{\link[aroma.affymetrix:extractExpressionSet.ChipEffectSet]{*extractExpressionSet}()}.
  For more information see \code{\link{AffymetrixCelSet}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
