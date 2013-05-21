%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  createExonByTranscriptCdf.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{createExonByTranscriptCdf.AffymetrixCdfFile}
\alias{createExonByTranscriptCdf.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.createExonByTranscriptCdf}
\alias{createExonByTranscriptCdf,AffymetrixCdfFile-method}

 \title{Creates an exon-by-transcript CDF}

 \description{
  Creates an exon-by-transcript CDF based on the probesets defined in an "exon-only" CDF
  and transcript-exon mapping of a NetAffx probeset annotation data file.
 }

 \usage{
\method{createExonByTranscriptCdf}{AffymetrixCdfFile}(cdf, csv, tags=c("*"), path=getPath(cdf),
  type=c("all", "core", "extended", "full", "main", "control", "cds"), subsetBy=NULL,
  within=NULL, ..., overwrite=FALSE, verbose=FALSE)
}

 \arguments{
  \item{cdf}{An \code{\link[aroma.affymetrix]{AffymetrixCdfFile}} specifying
     an "exon-only" CDF, which defines the exon-specific probesets
     that will go into the new CDF. For more details, see below.}
  \item{csv}{An \code{\link[aroma.affymetrix]{AffymetrixNetAffxCsvFile}}
     specifying the Affymetrix NetAffx CSV probeset annotation file
     that contains the transcript-exon mapping.}
  \item{tags}{Additional tags added to the filename of created CDF,
      i.e. <chiptype>,<tags>.cdf.}
  \item{path}{The output path where the custom CDF will be written.}
  \item{type}{A \code{\link[base]{character}} string specifying the type of CDF to be written.}
  \item{subsetBy}{An optional \code{\link[base]{character}} specifying the name of a column
     in the annotation file to subset against.  The column will be parsed
     as the data type of argument \code{within}.}
  \item{within}{A \code{\link[base]{vector}} of values accepted for the \code{subsetBy} column.}
  \item{...}{Additional arguments passed to \code{readDataFrame()} of
     \code{\link[aroma.affymetrix]{AffymetrixNetAffxCsvFile}}, e.g. \code{nrow}.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, an existing CDF is overwritten.}
  \item{verbose}{...}
 }

 \value{
   Returns an \code{\link[aroma.affymetrix]{AffymetrixCdfFile}}.
 }

 \section{Requirements for the "exon-only" CDF}{
   The template CDF - argument \code{cdf} - should be an "exon-only" CDF:
   each unit has one group/probeset, which is the exon.
   An example of this is the "unsupported" HuEx-1_0-st-v2.cdf
   from Affymetrix, which has 1,432,154 units.
   Such "exon-only" CDFs do not contain information about clustering
   exons/probesets into gene transcripts.
   The CDF may also contain a number of non-exon probesets corresponding
   to control probes, which can contain \emph{very} large numbers of
   probes per probeset. Such units are dropped/ignored by this method.
 }

 \section{Ordering of transcripts and exons within transcripts}{
   The transcripts (=units) will be ordered as they appear in the
   NetAffx annotation file.
   Within each transcript (=unit), the exons (=groups) are ordered
   lexicographically by their names.
   %% (Before Jan 28, 2008 (rev 3911) sorting was not done).
 }

 \section{Naming of transcripts and exons}{
   In the created CDF, each unit corresponds to one transcript cluster,
   and each group within a unit corresponds to the exons within
   that transcript cluster.  Thus, the unit names correspond to the
   transcript cluster names and the group names correspond to the
   exon names.

   The exon names are defined by unit names of the exon-only CDF,
   whereas the transcript names are defined by the
   \code{transcriptClusterId} column in the NetAffx annotation data file.
   These transcript and exon names are often "non-sense" integers.
   In order to map these to more genome-friendly names, use the various
   annotations available in the NetAffx annotation data file.
 }

 \examples{\dontrun{
 # The exon-only CDF
 cdf <- AffymetrixCdfFile$byChipType("HuEx-1_0-st-v2");

 # The NetAffx probeset annotation data file
 csv <- AffymetrixNetAffxCsvFile("HuEx-1_0-st-v2.na24.hg18.probeset.csv", path=getPath(cdf));

 # Create a CDF containing all core probesets:
 cdfT <- createExonByTranscriptCdf(cdf, csv=csv, tags=c("*,HB20110911"));
 print(cdfT);

 # Create CDF containing the core probesets with 3 or 4 probes:
 cdfT2 <- createExonByTranscriptCdf(cdf, csv=csv,
             tags=c("*,bySize=3-4,HB20110911"),
             subsetBy="probeCount", within=c("3", "4"));
 print(cdfT2);
 }}

 \author{
   Henrik Bengtsson adopted from \code{createTranscriptCDF()} written
   by Ken Simpson, Elizabeth Purdom and Mark Robinson.
 }


\keyword{internal}
\keyword{methods}
