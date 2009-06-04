%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getGenomeInformation.AffymetrixCdfFile}
\alias{getGenomeInformation.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.getGenomeInformation}
\alias{getGenomeInformation.AffymetrixCdfFile}
\alias{getGenomeInformation,AffymetrixCdfFile-method}


\title{Gets genome information for this chip type}

\description{
 Gets genome information for this chip type.
}

\usage{\method{getGenomeInformation}{AffymetrixCdfFile}(this, types=c("UGP", "dChip"), ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{types}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying what type of genome
    information sets to search for.}
  \item{...}{Not used.}
  \item{force}{If \code{\link[base:logical]{FALSE}}, cached information is retrieved, otherwise not.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link{GenomeInformation}} object.
}

\examples{\dontrun{
 
for (zzz in 0) {

# Define a CDF for a SNP chip
if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$fromChipType(".*Sty.*")
}
print(cdf)

# Get the genome information
gi <- getGenomeInformation(cdf)
print(gi)

# Get the units on chromosome X order position by default
# (First call will be slow be data is being cached)
units <- getUnitIndices(gi, chromosome="X")
nunits <- length(units)
cat(sprintf("Number of units on chromosome X: \%d\n", nunits))
str(units)

# Plot the SNP density for this chromosome
plotDensity(gi, chromosome="X");

} # for (zzz in 0)
rm(zzz)

}}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
