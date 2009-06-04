%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.extractAffyBatch.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{extractAffyBatch.AffymetrixCelSet}
\alias{extractAffyBatch.AffymetrixCelSet}
\alias{AffymetrixCelSet.extractAffyBatch}
\alias{extractAffyBatch.AffymetrixCelSet}
\alias{extractAffyBatch,AffymetrixCelSet-method}


\title{Extracts an in-memory AffyBatch object from the CEL set}

\description{
 Extracts an in-memory AffyBatch object from the CEL set.
 Note that any modifications done to the extract object will \emph{not}
 be reflected in the original CEL set.
}

\usage{\method{extractAffyBatch}{AffymetrixCelSet}(this, ..., verbose=FALSE)}

\arguments{
  \item{...}{Argument passed to \code{ReadAffy()}
    (\code{\link[affy]{read.affybatch}}).}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link[affy]{AffyBatch-class}} object.
}

\details{
 Since the \pkg{affy} package is making use of special CDF environment
 packages, this method will warn if the needed package is missing and
 explain that \pkg{affy} will later try to download and install it
 automatically.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
