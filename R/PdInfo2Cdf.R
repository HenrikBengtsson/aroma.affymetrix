
PdInfo2Cdf <- function( pdpkg, celfile, verbose=TRUE, overwrite=FALSE ) {

  #### This script has been written to generate a .cdf-file from an "pd.XXXX" package, 
  #### such as those build with pdInfoBuilder.
  #### The original was written by Samuel Wuest, modified by Mark Robinson 
  #### (around 12 Jan 2009) to be generic

  require("affxparser") || stop("Package not loaded: affxparser");
  require("pdInfoBuilder") || stop("Package not loaded: pdInfoBuilder");

  do.call("library",args=list(package=pdpkg))

  pmFeature2List <- function(u) {
    nr <- nrow(u)
    o <- order(u$atom)
    v <- 0:(nr-1)
    id <- u$fsetid[1]
    g <- list(list(x=u$x[o], y=u$y[o], pbase=rep("T",nr),tbase=rep("A",nr),
              atom=v,indexpos=v,groupdirection="sense",natoms=nr,ncellsperatom=1))
    names(g) <- id
    list(unittype=1,unitdirection=1,groups=g,natoms=nr,ncells=nr,ncellsperatom=1,unitnumber=id)
  }

  if(verbose)
    cat("Reading CEL file.\n")
  cel <- read.celfiles(filenames=celfile,pkgname=pdpkg)
  pd <-getPlatformDesign(cel)
  if(verbose)
    cat("Querying Platform Design database.\n")
  ff <- dbGetQuery(db(pd), "select * from pmfeature")

  celHead <- readCelHeader(celfile)
  nrows <- celHead$rows
  ncols <- celHead$rows


  # three 3 lines speed up the splitting ...
  if(verbose)
    cat("Creating list from query table.\n")
  ffs <- split(ff, substr(ff$fsetid,1,4))
  ffs <- unlist( lapply(ffs, FUN=function(u) split(u,u$fsetid)), recursive=FALSE)
  names(ffs) <- substr(names(ffs), 6,nchar(names(ffs)))

  nunits <- length(ffs)

  pdName <- gsub("\\.","",pdpkg)

  ## creating the cdf-header;
  newCdfHeader <- list(ncols=ncols, nrows=nrows, nunits=nunits, nqcunits=0, refseq="", 
                       chiptype=pdName, filename=paste(pdName,"cdf",sep="."), rows=nrows, 
                       cols=ncols, probesets= nunits, qcprobesets= 0, reference="")

  ### make the input-list for the writeCdf-function
  if(verbose)
    cat(paste("Creating CDF list for", nunits, "units.\n"))
  newCdfList <- lapply(ffs, pmFeature2List)

  ### writing the cdf-file (binary-file)
  writeCdf(newCdfHeader$filename, cdfheader= newCdfHeader, cdf=newCdfList, cdfqc=NULL, verbose=verbose, overwrite=overwrite)

}