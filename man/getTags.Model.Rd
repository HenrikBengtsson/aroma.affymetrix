%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Model.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getTags.Model}
\alias{getTags.Model}
\alias{Model.getTags}
\alias{getTags.Model}
\alias{getTags,Model-method}


\title{Gets the tags of the output data set}

\description{
 Gets the tags of the output data set, which consists of the tags of the input data set followed
 by an additional set of tags added by the model.
}

\usage{\method{getTags}{Model}(this, collapse=NULL, ...)}

\arguments{
 \item{collapse}{A \code{\link[base]{character}} string used to concatenate the tags.
    If \code{\link[base]{NULL}}, the tags are not concatenated.}
 \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{character}} \code{\link[base]{vector}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:setTags.Model]{*setTags}()}
  For more information see \code{\link{Model}}.
}
\keyword{internal}
\keyword{methods}
