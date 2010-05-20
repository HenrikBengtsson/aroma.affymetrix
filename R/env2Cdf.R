
Env2Cdf <- function( env, celfile, verbose=TRUE, overwrite=FALSE ) {

  #### This script has been written to generate a .cdf-file from an "XXXXcdf" package, 
  #### such as the Bioconductor 'metadata' packages.
  #### The original was written by Samuel Wuest, modified by Mark Robinson 
  #### (around 12 Jan 2009) to be generic

  require("affxparser") || stop("Package not loaded: affxparser");
  do.call("library",args=list(package=env))

  if(verbose)
    cat("Reading environment: ", env,".\n", sep="")
  ffs <- as.list(get(env))

  if(verbose)
    cat("Reading CEL file header.\n")
  celHead <- readCelHeader(celfile)
  nrows <- celHead$rows
  ncols <- celHead$rows

  env2List <- function(u) {
    pmi <- u[,"pm"]
    mmi <- u[,"mm"]
	
	x <- (u-1)%%nrows
	y <- (u-1)%/%nrows
    nr <- nrow(u)
	if( all(is.na(mmi)) ) {
	  # do PM only
	  # so far this is not implemented
	  #o <- which(!is.na(x))
      #g <- list(list(x=u$x[o], y=u$y[o], pbase=rep("T",nr),tbase=rep("A",nr),
      #        atom=v,indexpos=v,groupdirection="sense",natoms=nr,ncellsperatom=1))
      #names(g) <- id
      #list(unittype=1,unitdirection=1,groups=g,natoms=nr,ncells=nr,ncellsperatom=1,unitnumber=1)
	  return(NULL)
	} else {
	  # do PM and MM
      v <- 0:(nr-1)
	  atom <- rep(v,2)
	  o <- order(atom)
      g <- list(list(x=c(x)[o], y=c(y)[o], pbase=rep(c("A","T"),nr),tbase=rep(c("T","T"),nr),
              atom=atom[o],indexpos=atom[o],groupdirection="sense",natoms=nr,ncellsperatom=2))
      #names(g) <- id
      list(unittype=1,unitdirection=1,groups=g,natoms=nr,ncells=nr*2,ncellsperatom=2,unitnumber=1)
	}
  }

  nunits <- length(ffs)
  
  ### make the input-list for the writeCdf-function
  if(verbose)
    cat(paste("Creating CDF list for", nunits, "units.\n"))
  newCdfList <- lapply(ffs, env2List)
  #return(newCdfList)

  if(verbose)
    cat("Adding group names to CDF list.\n")
  nm <- names(newCdfList)
  for(i in 1:nunits) {
    names( newCdfList[[i]]$groups ) <- nm[i]
    newCdfList[[i]]$unitnumber <- i
  }

  pdName <- gsub("cdf","",env)
  ## creating the cdf-header;
  newCdfHeader <- list(ncols=ncols, nrows=nrows, nunits=nunits, nqcunits=0, refseq="", 
                       chiptype=pdName, filename=paste(pdName,"cdf",sep="."), rows=nrows, 
                       cols=ncols, probesets= nunits, qcprobesets= 0, reference="")


  ### writing the cdf-file (binary-file)
  writeCdf(newCdfHeader$filename, cdfheader= newCdfHeader, cdf=newCdfList, cdfqc=NULL, verbose=verbose, overwrite=overwrite)

}