setMethodS3("createBlankSet", "AffymetrixCelSet", function(static, name, tags=NULL, sampleNames, celTemplate, path="cel", overwrite=FALSE, ..., verbose=FALSE) {
  # Argument 'celTemplate':
  celTemplate <- Arguments$getInstanceOf(celTemplate, "AffymetrixCelFile");

  # Argument 'name':
  name <- Arguments$getCharacter(name, length=c(1,1));

  # Argument 'sampleNames':
  sampleNames <- Arguments$getCharacters(sampleNames);

  # Argument 'tags':
  tags <- Arguments$getCharacters(tags);

  # Argument 'path':
  path <- Arguments$getWritablePath(path);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the CDF
  cdf <- getCdf(celTemplate);
  chipType <- getChipType(cdf);

  # The data-set path
  fullname <- paste(c(name, tags), collapse=",");
  path <- filePath(path, fullname, chipType, expandLinks="any");
  mkdirs(path);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Create output pathnames and check overwriting rules
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  pathnames <- c();
  nbrOfFiles <- length(sampleNames);
  for (kk in seq(length=nbrOfFiles)) {
    filename <- sprintf("%s.CEL", sampleNames[kk]);
    pathname <- Arguments$getWritablePathname(filename, path=path,
                                                  mustNotExist=!overwrite);
    pathnames[kk] <- pathname;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Create CEL files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdfHeader <- getHeader(cdf);
  celPathname <- getPathname(celTemplate);
  celHeader0 <- readCelHeader(celPathname);
 
  for (kk in seq(length=nbrOfFiles)) {
    pathname <- pathnames[kk];
    fullname <- basename(pathname);
    fullname <- gsub("[.](c|C)(e|E)(l|L)$", "", fullname);
    verbose && enter(verbose, sprintf("Array #%d of %d", kk, nbrOfFiles));

    verbose && cat(verbose, "Full name: ", fullname);
    verbose && cat(verbose, "Pathname: ", pathname);

    # Build a valid CEL header
    celHeader <- cdfHeaderToCelHeader(cdfHeader, sampleName=fullname);
  
    # Add some extra information about what the CEL file is for
    params <- c(Descripion="This CEL file contains probe data saved by the aroma.affymetrix package.");
    parameters <- gsub(" ", "_", params);
    names(parameters) <- names(params);
    parameters <- paste(names(parameters), parameters, sep=":");
    parameters <- paste(parameters, collapse=";");
    parameters <- paste(celHeader$parameters, parameters, "", sep=";");
    parameters <- gsub(";;", ";", parameters);
    parameters <- gsub(";$", "", parameters);
    celHeader$parameters <- parameters;
    rm(params, parameters);

    # Create an empty CEL file
    verbose && enter(verbose, "Creating empty CEL file");
    createCel(pathname, header=celHeader, overwrite=overwrite, ..., 
                                                     verbose=less(verbose));
    rm(celHeader);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (kk ...)

  # Finally, return the AffymetrixCelSet object
  byPath(static, path);
}, static=TRUE, private=TRUE)



############################################################################
# HISTORY:
# 2007-08-09
# o AffymetrixCelSet$createBlankSet() now creates CEL files with upper-case
#   filename extension "*.CEL", not "*.cel".  The reason for this is that
#   some software don't recognize lower case filename extensions :(
#   Note: The above modification is not safe in the sense that if CEL files
#   in lower-case already exists, the above method will not detect those
#   on case-sensitive file systems (e.g. Unix), and create a new set.
#   However, that should be fine because this functions has in practice not
#   been used by anyone.
# 2006-12-06
# o Created.
############################################################################
