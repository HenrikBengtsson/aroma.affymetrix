%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaChipTypeAnnotationFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fromFile.AromaChipTypeAnnotationFile}
\alias{fromFile.AromaChipTypeAnnotationFile}
\alias{AromaChipTypeAnnotationFile.fromFile}
\alias{fromFile.AromaChipTypeAnnotationFile}
\alias{fromFile,AromaChipTypeAnnotationFile-method}


\title{Sets up an AromaChipTypeAnnotationFile}

\description{
 Sets up an AromaChipTypeAnnotationFile.
}

\synopsis{fromFile.AromaChipTypeAnnotationFile(static, filename, path=NULL, ...)}
\usage{AromaChipTypeAnnotationFile$fromFile(filename, path=NULL, ...)}

\arguments{
 \item{filename}{The filename of to the file.}
 \item{path}{The path to the file.}
 \item{...}{Not used.}
}

\value{
 Returns an instance of \code{\link{AromaChipTypeAnnotationFile}} or its subclasses.
 If the file is not found or if it is of the wrong file format, an
 error is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:byChipType.AromaChipTypeAnnotationFile]{*byChipType}()}.
  For more information see \code{\link{AromaChipTypeAnnotationFile}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
