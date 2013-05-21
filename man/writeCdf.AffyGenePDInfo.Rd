%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffyGenePDInfo.writeCdf.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{writeCdf.AffyGenePDInfo}
\alias{writeCdf.AffyGenePDInfo}
\alias{AffyGenePDInfo.writeCdf}
\alias{writeCdf,AffyGenePDInfo-method}
 \alias{writeCdf.PDInfoList}
 \alias{writeCdf.DBPDInfo}

 \title{Generates an Affymetrix CDF file from a Platform Design (PD) package}

 \description{
   Generates an Affymetrix CDF file from a Platform Design (PD) package.
   Platform Design (PD) packages are also known as "pdInfo" packages.
 }

 \usage{
\method{writeCdf}{AffyGenePDInfo}(this, tags=c("*"), unitsBy=c("transcript", "exon"), namesBy=c("fsetid", "id"),
  path=NULL, overwrite=FALSE, verbose=TRUE, ...)
}

 \arguments{
  \item{tags}{An optional \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be added to the CDF
    filename.}
  \item{unitsBy}{A \code{\link[base]{character}} string specifying how to group units.}
  \item{path}{The output path where the CDF file is written.
    If \code{\link[base]{NULL}} (default), then it is written to the corresponding
    \code{annotationData/chipTypes/<chipType>/} directory.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, an existing CDF is overwritten, otherwise
    an exception is thrown.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
  \item{...}{Not used.}
 }

 \value{
   Returns (invisibly) the pathname to CDF written.
 }

 \details{
   The formal chip type of the CDF is inferred from the AffyGenePDInfo package.
   The filename of the CDF is constructed from the chip type and any optional
   tags.
   To minimize the risk for a corrupt CDF file, the creation of the file
   is atomic by first writing to a temporary file which is then renamed.
 }

 \section{Limitations}{
   The information available in the PD package is limited and does
   not contain all information needed to populate a CDF file.
   In order to workaround these limitations, certain CDF entries
   are set to predefined/hardwired values.
   The 'pbase' and 'tbase' entries of the generated CDF file is
   hardwired to "T" and "A", respectively.
   Likewise, the 'groupdirection' entry is hardwired to "sense".
 }

 \author{
   Henrik Bengtsson and Guido Hooiveld adopted from \code{pdInfo2Cdf()}
   written by Samuel Wuest and Mark Robinson.
 }


\keyword{internal}
\keyword{methods}
