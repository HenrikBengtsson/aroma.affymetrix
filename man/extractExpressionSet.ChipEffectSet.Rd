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
\alias{extractExpressionSet.ChipEffectSet}
\alias{extractExpressionSet,ChipEffectSet-method}


\title{Extracts an in-memory ExpressionSet object}

\description{
 Extracts an in-memory ExpressionSet object from a \code{\link{ChipEffectSet}} object.
 Note that any modifications done to the extract object will \emph{not}
 be reflected in the original chip-effect set.
}

\usage{\method{extractExpressionSet}{ChipEffectSet}(this, ..., logBase=2, verbose=FALSE)}

\arguments{
  \item{...}{Optional argument passed to \code{\link[aroma.affymetrix:extractMatrix.ChipEffectSet]{*extractMatrix}()}.}
  \item{logBase}{An \code{\link[base]{integer}} specifying the base to use when
   log-transforming the signals.  If \code{\link[base]{NULL}}, the signals are not
   transformed, but kept as is.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link[Biobase]{ExpressionSet-class}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ChipEffectSet}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
