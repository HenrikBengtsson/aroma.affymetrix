PdInfo2Cdf <- function(pdpkg, celfile, verbose=TRUE, overwrite=FALSE) {
  #### This script has been written to generate a .cdf-file from an 
  #### "pd.XXXX" package, such as those build with pdInfoBuilder.
  #### The original was written by Samuel Wuest, modified by Mark Robinson 
  #### (around 12 Jan 2009) to be generic

  require("affxparser") || throw("Package not loaded: affxparser");
  require("pdInfoBuilder") || throw("Package not loaded: pdInfoBuilder");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pmFeature2List <- function(u) {
    nr <- nrow(u);
    o <- order(u$atom);
    v <- 0:(nr-1);
    id <- u$fsetid[1];
    g <- list(list(x=u$x[o], y=u$y[o], pbase=rep("T", nr), tbase=rep("A", nr),
              atom=v, indexpos=v, groupdirection="sense", natoms=nr,
              ncellsperatom=1));
    names(g) <- id;
    list(unittype=1, unitdirection=1, groups=g, natoms=nr, ncells=nr,
         ncellsperatom=1, unitnumber=id);
  } # pmFeature2List()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'celfile':
  celfile <- Arguments$getReadablePathname(celfile, mustExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);


  verbose && enter(verbose, "PdInfo2Cdf()");

  filename <- sprintf("%s.cdf", pdName);
  filename <- Arguments$getWritablePathname(filename, mustNotExist=!overwrite);

  # Load the required PD package.
  require(pdpkg, character.only=TRUE) || 
                 throw("Platform Design (PD) package not loaded: ", pdpkg);

  verbose && cat(verbose, "Reading CEL file.");

  verbose && enter(verbose, "Reading CEL file");
  cel <- read.celfiles(filenames=celfile, pkgname=pdpkg);
  verbose && exit(verbose);

  pd <- getPlatformDesign(cel);
  rm(cel);  # Not needed anymore

  verbose && enter(verbose, "Querying Platform Design database");
  ff <- dbGetQuery(db(pd), "select * from pmfeature");
  rm(pd);  # Not needed anymore
  verbose && exit(verbose);

  celHead <- readCelHeader(celfile);
  nrows <- celHead$rows;
  ncols <- celHead$rows;
  rm(celHead);  # Not needed anymore

  verbose && enter(verbose, "Creating list from query table");
  # three 3 lines speed up the splitting ...
  ffs <- split(ff, substr(ff$fsetid, 1,4));
  ffs <- lapply(ffs, FUN=function(u) split(u, u$fsetid));
  ffs <- unlist(ffs, recursive=FALSE);
  names(ffs) <- substr(names(ffs), 6, nchar(names(ffs)));
  verbose && exit(verbose);

  nunits <- length(ffs);
  pdName <- gsub("\\.", "", pdpkg);

  ## creating the CDF header;
  newCdfHeader <- list(ncols=ncols, nrows=nrows, nunits=nunits, 
                       nqcunits=0, refseq="", chiptype=pdName, 
                       filename=filename, rows=nrows, 
                       cols=ncols, probesets= nunits, 
                       qcprobesets=0, reference="");

  verbose && enter(verbose, "Creating CDF list structure");
  ### make the input-list for writeCdf()
  verbose && cat(verbose, "Creating CDF list for ", nunits, " units");
  newCdfList <- lapply(ffs, FUN=pmFeature2List);
  rm(ffs);  # Not needed anymore
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing CDF file");
  ### writing the CDF file (binary-file)
  res <- writeCdf(newCdfHeader$filename, cdfheader=newCdfHeader, 
          cdf=newCdfList, cdfqc=NULL, verbose=verbose, overwrite=overwrite);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(res);
} # PdInfo2Cdf()

############################################################################
# HISTORY:
# 2009-10-16 [HB]
# o MEMORY OPTIMIZATION: Cleaning out non needed objects sooner.
# o Added verbose a'la R.utils.
# o Added some validation of arguments.
# o Tidied up the code structure.
# 2009-01-13 [MR]
# o Added.
# 2008-??-?? [SW]
# o Created by Samuel Wuest.
############################################################################
