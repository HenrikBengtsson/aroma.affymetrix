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

\usage{
\method{bpmapCluster2Cdf}{default}(pathname, chipType, tags=NULL, rows, cols, maxProbeDistance=3000L,
  minNbrOfProbes=30L, groupName=gsub("_.*", "", chipType), field="fullname",
  stringRemove=sprintf("\%s:.*;", groupName), ..., flavor=c("v2", "v1"), path="*",
  verbose=-10)
}

\arguments{
 \item{pathname}{The pathname to the BPMAP file.}
 \item{chipType, tags}{The chip type and optional tags of the CDF to
   be written.}
 \item{rows, cols}{Two positive \code{\link[base]{integer}}s specifying the probe dimension
   of the chip.  It is important to get this correct.  They can be
   inferred from the CEL header of a CEL file for this chip,
   cf. \code{\link[affxparser]{readCelHeader}}.}
 \item{maxProbeDistance}{A positive \code{\link[base]{integer}} specifying the maximum
   genomic distance (in basepairs) allowed between two probes in order
   to "cluster" those two probes into the same CDF units.  Whenever the
   distance is greater, the two probes end up in two different CDF units.}
 \item{minNbrOfProbes}{A positive \code{\link[base]{integer}} specifying the minimum number
   of probes required in a CDF unit.  If fewer, those probes are dropped.}
 \item{groupName}{A \code{\link[base]{character}} string specifying which BPMAP sequences
    to keep.  Sequence with this group name is kept, all others are
    excluded.}
 \item{field}{A \code{\link[base]{character}} string.}
 \item{stringRemove}{An (optional) regular expression.}
 \item{...}{Optional arguments passed to \code{\link[affxparser]{readBpmap}}.}
 \item{flavor}{Specifying which version of BPMAP-to-CDF generator
    to use. The default is always to use the most recent one, which
    is also the recommended one.  Previous versions are kept only for
    backward compatibility (and may be dropped at anytime).}
 \item{path}{The directory where the CDF file will be written.
    If \code{"*"} (default), it will be written to
    \code{annotationData/chipTypes/<chipType>/}.}
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
  Flavor \code{"v2"} replaced \code{"v1"} as aroma.affymetrix v2.5.4
  (June 21, 2012). For details, see \code{news(package="aroma.affymetrix")}.
}

\author{
  Henrik Bengtsson adopted from Mark Robinson standalone/online version
  as of July 11, 2011.
}


\keyword{internal}
