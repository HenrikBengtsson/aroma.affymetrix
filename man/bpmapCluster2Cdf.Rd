%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  bpmapCluster2Cdf.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{bpmapCluster2Cdf}
\alias{bpmapCluster2Cdf.default}
\alias{bpmapCluster2Cdf}

 \title{Creates a CDF from tiling-array BPMAP file}

 \description{
   Creates a CDF from tiling-array BPMAP file.\cr

   \emph{
    NOTE: This method applies only to Affymetrix tiling arrays!
    Furthermore, it is likely to be more useful for promoter tiling arrays
    and less so for whole-genome tiling arrays.
   }
 }

 \usage{\method{bpmapCluster2Cdf}{default}(pathname, chipType, tags=NULL, maxProbeDistance=3000, minNbrOfProbes=30, ..., rows=NULL, cols=NULL, groupName=gsub("_.*", "", chipType), field="fullname", stringRemove=sprintf("\%s:.*;", groupName), verbose=-10)}

 \arguments{
  \item{pathname}{The pathname to the BPMAP file.}
  \item{chipType, tags}{The chip type and optional tags of the CDF to
    be written.}
  \item{maxProbeDistance}{A positive \code{\link[base]{integer}} specifying the maximum
    genomic distance (in basepairs) allowed between two probes in order
    to "cluster" those two probes into the same CDF units.  Whenever the
    distance is greater, the two probes end up in two different CDF units.}
  \item{minNbrOfProbes}{A positive \code{\link[base]{integer}} specifying the minimum number
    of probes required in a CDF unit.  If fewer, those probes are dropped.}
  \item{...}{Not used.}
  \item{rows, cols}{Two (optional) positive \code{\link[base]{integer}}s.
     If \code{\link[base]{NULL}}, optimal values are inferred auotmatically.}
  \item{groupName}{A \code{\link[base]{character}} string.}
  \item{field}{A \code{\link[base]{character}} string.}
  \item{stringRemove}{An (optional) regular expression.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
  Returns (invisibly) a the pathname of the created CDF file.
  The created CDF is written to the current directory.
 }

 \details{
   This method applies only to Affymetrix tiling arrays.  It is likely
   to be useful for promoter tiling arrays and less so for whole-genome
   tiling arrays.
 }

 \author{
   Henrik Bengtsson adopted from Mark Robinson standalone/online version
   as of July 11, 2011.
 }


\keyword{internal}
