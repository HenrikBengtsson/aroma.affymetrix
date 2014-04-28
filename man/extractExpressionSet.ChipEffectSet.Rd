%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectSet.extractExpressionSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{extractExpressionSet.ChipEffectSet}
\alias{extractExpressionSet.ChipEffectSet}
\alias{ChipEffectSet.extractExpressionSet}
\alias{extractExpressionSet,ChipEffectSet-method}
\alias{extractExpressionSet}

\title{Extracts an in-memory ExpressionSet object}

\description{
 Extracts an in-memory ExpressionSet object from a \code{\link{ChipEffectSet}} object.
 Note that any modifications done to the extracted object will \emph{not}
 be reflected in the original chip-effect set.
}

\usage{
\method{extractExpressionSet}{ChipEffectSet}(this, ..., logBase=2, orderUnitsBy=c("asis", "lexicographic"),
  annotationPkg=NULL, verbose=FALSE)
}

\arguments{
  \item{...}{Optional argument passed to \code{\link[aroma.affymetrix:extractMatrix.ChipEffectSet]{*extractMatrix}()}.}
  \item{logBase}{An \code{\link[base]{integer}} specifying the base to use when
   log-transforming the signals.  If \code{\link[base]{NULL}}, the signals are not
   transformed, but kept as is.}
  \item{annotationPkg}{(optional) A \code{\link[base]{character}} string specifies a
   Bioconductor (ChipDb, CDF environment or PDInfo) annotation package,
   which then sets the 'annotation' slot of the returned object.
   If \code{"ChipDb"}, \code{"cdf"} or \code{"PDInfo"} th corresponding
   ChipDB, CDF environment or PDInfo package is inferred from
   the CDF's chip type.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link[Biobase]{ExpressionSet-class}} object.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{ChipEffectSet}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
