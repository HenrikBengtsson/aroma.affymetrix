%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnagCfhFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromFile.CnagCfhFile}
\alias{fromFile.CnagCfhFile}
\alias{CnagCfhFile.fromFile}
\alias{fromFile.CnagCfhFile}
\alias{fromFile,CnagCfhFile-method}


\title{Defines an CnagCfhFile object from a CFH file}

\description{
 Defines an CnagCfhFile object from a CFH file.
}

\synopsis{fromFile.CnagCfhFile(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE)}
\usage{CnagCfhFile$fromFile(filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE)}

\arguments{
 \item{filename}{The filename of to the file.}
 \item{path}{The path to the file.}
 \item{...}{Not used.}
}

\value{
 Returns an \code{\link{CnagCfhFile}} object.
 If the file is not found or if it is of the wrong file format, an
 error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{CnagCfhFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
