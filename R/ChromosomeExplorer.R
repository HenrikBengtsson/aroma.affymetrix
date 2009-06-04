###########################################################################/**
# @RdocClass ChromosomeExplorer
#
# @title "The ChromosomeExplorer class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{model}{A @see "CopyNumberChromosomalModel" object.}
#   \item{zooms}{An positive @integer @vector specifying for which zoom
#    levels the graphics should be generated.}
#   \item{...}{Not used.}
#   \item{version}{The version of the Explorer HTML/Javascript generated/used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Generating PNG images}{
#   In order to get better looking graphs, but also to be able to generate
#   bitmap images on systems without direct bitmap support, which is the case
#   when running R in batch mode or on Unix without X11 support, images are
#   created using the @see "R.utils::png2" device (a wrapper for 
#   \code{bitmap()} immitating \code{png()}).  The \code{png()} is only
#   used if \code{png2()}, which requires Ghostscript, does not.
#   Note, when images are created using \code{png2()}, the images does
#   not appear immediately, although the function call is completed,
#   so be patient.
# }
#
# @author
# 
# \seealso{
#  @see "CopyNumberChromosomalModel".
# }
#*/###########################################################################
setConstructorS3("ChromosomeExplorer", function(model=NULL, zooms=2^(0:7), ..., version=c("3")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  if (!is.null(model)) {
    if (!inherits(model, "CopyNumberChromosomalModel")) {
      throw("Argument 'model' is not a 'CopyNumberChromosomalModel': ", 
                                                            class(model)[1]);
    }
  }

  # Argument 'zooms':
  if (is.null(zooms)) {
    zooms <- 2^(0:7);
  } else {
    zooms <- Arguments$getDoubles(zooms, range=c(0, Inf));
  }
  
  # Argument 'version':
  version <- match.arg(version);


  extend(Explorer(...), "ChromosomeExplorer",
    .model = model,
    .arrays = NULL,
    .plotCytoband = TRUE,
    .zooms = zooms,
    .version = version
  )
})


setMethodS3("as.character", "ChromosomeExplorer", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Number of arrays:", nbrOfArrays(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("setCytoband", "ChromosomeExplorer", function(this, status=TRUE, ...) {
  # Argument 'status':
  status <- Arguments$getLogical(status);

  this$.plotCytoband <- status;
})


###########################################################################/**
# @RdocMethod getModel
#
# @title "Gets the model"
#
# \description{
#  @get "title" for which the explorer is displaying it results.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "CopyNumberChromosomalModel".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getModel", "ChromosomeExplorer", function(this, ...) {
  this$.model;
})


setMethodS3("getNames", "ChromosomeExplorer", function(this, ...) {
  getArrays(this, ...);
})

setMethodS3("getArraysOfInput", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getNames(model, ...);
}, protected=TRUE)


setMethodS3("getFullNames", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  arrays <- getArrays(this);
  idx <- match(arrays, getNames(model));
  fullnames <- getFullNames(model, arrays=idx, ...);
  fullnames;
})


###########################################################################/**
# @RdocMethod setArrays
#
# @title "Sets the arrays"
#
# \description{
#  @get "title" to be processed by the explorer.
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @character (or @integer) @vector of arrays used in the
#      model.  If @NULL, all arrays of the model are considered.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setArrays", "ChromosomeExplorer", function(this, arrays=NULL, ...) {
  # Argument 'arrays':
  if (!is.null(arrays)) {
    setTuple <- getSetTuple(this);
    arrays <- indexOfArrays(setTuple, arrays);
    arrays <- getArrays(setTuple)[arrays];
  }

  this$.arrays <- arrays;
})

setMethodS3("getSetTuple", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getSetTuple(model);
}, protected=TRUE)


setMethodS3("getNameOfInput", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getName(model, ...);
}, protected=TRUE)


setMethodS3("getTagsOfInput", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getTags(model);
}, protected=TRUE)


setMethodS3("getPath", "ChromosomeExplorer", function(this, ...) {
  mainPath <- getMainPath(this);

  # Chip type    
  model <- getModel(this);
  chipType <- getChipType(model);

  # Image set
  set <- getSetTag(model);

  # The full path
  path <- filePath(mainPath, chipType, set, expandLinks="any");

  # Create path?
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})


###########################################################################/**
# @RdocMethod getChromosomes
#
# @title "Gets the chromosomes available"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChromosomes", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getChromosomes(model);
})


setMethodS3("getZooms", "ChromosomeExplorer", function(this, ...) {
  zooms <- this$.zooms;
  if (is.null(zooms))
    zooms <- 2^(0:7);
  zooms <- as.integer(zooms);
  zooms;
})


setMethodS3("setZooms", "ChromosomeExplorer", function(this, zooms=NULL, ...) {
  # Argument 'zooms':
  if (is.null(zooms)) {
    zooms <- 2^(0:7);
  } else {
    zooms <- Arguments$getIntegers(zooms, range=c(1, Inf));
  }
  zooms <- unique(zooms);
  zooms <- as.integer(zooms);
  this$.zooms <- zooms;
  invisible(this);
})




###########################################################################/**
# @RdocMethod updateSamplesFile
#
# @title "Updates the samples.js file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns (invisibly) the pathname to the samples.js file.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("updateSamplesFile", "ChromosomeExplorer", function(this, ..., verbose=FALSE) {
  require("R.rsp") || throw("Package not loaded: R.rsp");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Output version
  version <- this$.version;
  verbose && cat(verbose, "Explorer output version: ", version);


  path <- getPath(this);
  parentPath <- getParent(path);
  parent2Path <- getParent(parentPath);
  parent3Path <- getParent(parent2Path);

  srcPath <- getTemplatePath(this);
  if (version == "3") {
    pathname <- filePath(srcPath, "rsp", "ChromosomeExplorer3", "ChromosomeExplorer.onLoad.js.rsp");
  } else if (version == "4") {
    pathname <- filePath(srcPath, "rsp", "ChromosomeExplorer4", "ChromosomeExplorer4.onLoad.js.rsp");
  } else if (version == "5") {
    pathname <- filePath(srcPath, "rsp", "ChromosomeExplorer5", "ChromosomeExplorer5.onLoad.js.rsp");
  }
  verbose && enter(verbose, "Compiling ", basename(pathname));
  verbose && cat(verbose, "Source: ", pathname);
  outFile <- gsub("[.]rsp$", "", basename(pathname));
  outPath <- parent2Path;
  outPath <- Arguments$getWritablePath(outPath); 
  verbose && cat(verbose, "Output path: ", outPath);

  verbose && enter(verbose, "Scanning directories for available chip types");
  # Find all directories matching these
  dirs <- list.files(path=parent2Path, full.names=TRUE);
  dirs <- dirs[sapply(dirs, FUN=isDirectory)];

  # Get possible chip types
  model <- getModel(this);
  chipTypes <- c(getChipType(model), getChipTypes(model));
  chipTypes <- intersect(chipTypes, basename(dirs));
  verbose && cat(verbose, "Detected chip types: ", 
                                           paste(chipTypes, collapse=", "));
  verbose && exit(verbose);

  # Get available zooms
  verbose && enter(verbose, "Scanning image files for available zooms");
  pattern <- ".*,x([0-9][0-9]*)[.]png$";
  zooms <- list.files(path=path, pattern=pattern);
  zooms <- gsub(pattern, "\\1", zooms);
  zooms <- gsub("^0*", "", zooms);
  if (length(zooms) == 0) {
    # Default zooms
    zooms <- getZooms(this);
  }
  zooms <- unique(zooms);
  zooms <- as.integer(zooms);
  zooms <- sort(zooms);
  verbose && cat(verbose, "Detected (or default) zooms: ", paste(zooms, collapse=", "));
  verbose && exit(verbose);

  # Get available subdirectories
  verbose && enter(verbose, "Scanning directory for subdirectories");
  subdirs <- list.files(path=parentPath, full.names=TRUE);
  # Keep only directories
  subdirs <- subdirs[sapply(subdirs, FUN=isDirectory)];

  isChrLayer <- (regexpr(",chrLayer", subdirs) != -1);
  isSampleLayer <- (regexpr(",sampleLayer", subdirs) != -1);
  isLayer <- (isChrLayer | isSampleLayer);
  hasLayers <- any(isLayer);

  # Remove anything that is a layer
  sets <- basename(subdirs)[!isLayer];

  if (length(sets) == 0) {
    sets <- c("glad", "cbs");
    warning("No 'set' directories detected. Using defauls: ", paste(sets, collapse=", "));
  } else {
    sets <- basename(sets);
  }
  sets <- sort(unique(sets));
  verbose && cat(verbose, "Detected (or default) sets: ", paste(sets, collapse=", "));
  verbose && exit(verbose);

  if (hasLayers) {
    sets <- c("-LAYERS-", sets);
  }

  chrLayers <- gsub(",chrLayer", "", basename(subdirs)[isChrLayer]);
  sampleLayers <- gsub(",sampleLayer", "", basename(subdirs)[isSampleLayer]);

  # Compile RSP file
  verbose && enter(verbose, "Compiling RSP");
  env <- new.env();
  env$chipTypes <- chipTypes;
  env$samples <- getFullNames(this, ...);
  env$sampleLabels <- getNames(this);
  env$zooms <- zooms;
  env$sets <- sets;
  env$chrLayers <- chrLayers;
  env$sampleLayers <- sampleLayers;

  if (getParallelSafe(this)) {
    tryCatch({
      pathname <- rspToHtml(pathname, path=NULL, 
                            outFile=outFile, outPath=outPath, 
                            overwrite=TRUE, envir=env);
    }, error = function(ex) {})
  } else {
    pathname <- rspToHtml(pathname, path=NULL, 
                          outFile=outFile, outPath=outPath, 
                          overwrite=TRUE, envir=env);
  }

  verbose && exit(verbose);


  verbose && exit(verbose);
  
  invisible(pathname);
}, private=TRUE)



setMethodS3("addIndexFile", "ChromosomeExplorer", function(this, filename=NULL, ...) {
  if (is.null(filename)) {
    version <- this$.version;
    if (version == "3") {
      filename <- "ChromosomeExplorer.html";
    } else if (version == "4") {
      filename <- "ChromosomeExplorer4.html";
    } else if (version == "5") {
      filename <- "ChromosomeExplorer5.html";
    }
  }
  NextMethod("addIndexFile", this, filename=filename, ...);
}, protected=TRUE)




setMethodS3("setup", "ChromosomeExplorer", function(this, ..., force=FALSE) {
  if (getParallelSafe(this)) {
    tryCatch({
      # Setup includes/?
      addIncludes(this, ..., force=force);
  
      # Setup main HTML file
      addIndexFile(this, ..., force=force);
  
      # Setup samples.js?
      updateSamplesFile(this, ...);
    }, error = function(ex) {})
  } else {
    # Setup includes/?
    addIncludes(this, ..., force=force);

    # Setup main HTML file
    addIndexFile(this, ..., force=force);

    # Setup samples.js?
    updateSamplesFile(this, ...);
  }
}, private=TRUE)




setMethodS3("writeAxesLayers", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);

  path <- getPath(this);
  path <- filePath(getParent(path), "axes,chrLayer");


  if (getParallelSafe(this)) {
    tryCatch({
      path <- Arguments$getWritablePath(path);
      plotAxesLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);
    }, error = function(ex) {})
  } else {
    path <- Arguments$getWritablePath(path);
    plotAxesLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);
  }

  invisible(path);
}, private=TRUE)


setMethodS3("writeGridHorizontalLayers", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);

  path <- getPath(this);
  path <- filePath(getParent(path), "gridH,chrLayer");


  if (getParallelSafe(this)) {
    tryCatch({
      path <- Arguments$getWritablePath(path);
      plotGridHorizontalLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);
    }, error = function(ex) {})
  } else {
    path <- Arguments$getWritablePath(path);
    plotGridHorizontalLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);
  }

  invisible(path);
}, private=TRUE)



setMethodS3("writeCytobandLayers", "ChromosomeExplorer", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model
  model <- getModel(this);

  path <- getPath(this);
  path <- filePath(getParent(path), "cytoband,chrLayer");


  if (getParallelSafe(this)) {
    tryCatch({
      path <- Arguments$getWritablePath(path);
      plotCytobandLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);
    }, error = function(ex) {})
  } else {
    path <- Arguments$getWritablePath(path);
    plotCytobandLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);
  }

  invisible(path);
}, private=TRUE)


setMethodS3("getSampleLayerName", "Explorer", function(this, name, class="sampleLayer", ...) {
  layer <- c(getSampleLayerPrefix(this), name, class);
  layer <- paste(layer, collapse=",");
  layer;
}, private=TRUE)


setMethodS3("writeRawCopyNumberLayers", "ChromosomeExplorer", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model
  model <- getModel(this);

  path <- getPath(this);
  layer <- getSampleLayerName(this, "rawCNs");
  path <- filePath(getParent(path), layer);

  path <- Arguments$getWritablePath(path);
  plotRawCopyNumbers(model, path=path, imageFormat="png", transparent=TRUE, ...);

  invisible(path);
}, private=TRUE)


setMethodS3("writeCopyNumberRegionLayers", "ChromosomeExplorer", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model
  model <- getModel(this);

  # Not supported?
  if (!inherits(model, "CopyNumberSegmentationModel"))
    return(NULL);
    
  path <- getPath(this);
  tags <- getAsteriskTags(model, collapse=",");
  path <- filePath(getParent(path), sprintf("%s,sampleLayer", tags));

  path <- Arguments$getWritablePath(path);
  plotCopyNumberRegionLayers(model, path=path, imageFormat="png", transparent=TRUE, ...);

  invisible(path);
}, private=TRUE)




setMethodS3("writeGraphs", "ChromosomeExplorer", function(x, arrays=NULL, ...) {
  # To please R CMD check.
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (is.null(arrays))
    arrays <- getNames(this);


  # Get the model
  model <- getModel(this);

  path <- getPath(this);
  path <- Arguments$getWritablePath(path); 

  plotband <- this$.plotCytoband;  # Plot cytoband?
  plot(model, path=path, imageFormat="png", plotband=plotband, arrays=arrays, ...);

  invisible(path);
}, private=TRUE)



setMethodS3("writeRegions", "ChromosomeExplorer", function(this, arrays=NULL, nbrOfSnps=c(3,Inf), smoothing=c(-Inf,-0.15, +0.15,+Inf), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (is.null(arrays)) 
    arrays <- getNames(this);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  model <- getModel(this);

  # Not supported?
  if (!inherits(model, "CopyNumberSegmentationModel"))
    return(NULL);

  verbose && enter(verbose, "Writing CN regions");

  path <- getPath(this);
  path <- Arguments$getWritablePath(path); 

  dest <- filePath(path, "regions.xls");
  dest <- Arguments$getWritablePathname(dest); 


  # Extract and write regions
  if (getParallelSafe(this)) {
    tryCatch({
      pathname <- writeRegions(model, arrays=arrays, nbrOfSnps=nbrOfSnps, smoothing=smoothing, ..., skip=FALSE, verbose=less(verbose));
      res <- copyFile(pathname, dest, overwrite=TRUE);
      if (!res)
        dest <- NULL;
    }, error = function(ex) {})
  } else {
    pathname <- writeRegions(model, arrays=arrays, nbrOfSnps=nbrOfSnps, smoothing=smoothing, ..., skip=FALSE, verbose=less(verbose));
    res <- copyFile(pathname, dest, overwrite=TRUE);
    if (!res)
      dest <- NULL;
  }


  verbose && exit(verbose);

  invisible(dest);
}, private=TRUE)



###########################################################################/**
# @RdocMethod process
#
# @title "Generates image files, scripts and dynamic pages for the explorer"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @vector of arrays specifying which arrays to
#    be considered.  If @NULL, all are processed.}
#   \item{chromosome}{A @vector of chromosomes specifying which
#     chromosomes to be considered.  If @NULL, all are processed.}
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "ChromosomeExplorer", function(this, arrays=NULL, chromosomes=NULL, ..., zooms=getZooms(this), layers=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (is.null(arrays))
    arrays <- getNames(this);

  # Argument 'chromosomes':
  allChromosomes <- getChromosomes(this);
  if (is.null(chromosomes)) {
    chromosomes <- allChromosomes;
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  
  verbose && enter(verbose, "Generating ChromosomeExplorer report");

  # Setup HTML, CSS, Javascript files first
  if (getParallelSafe(this)) {
    tryCatch({
      setup(this, ..., verbose=less(verbose));
    }, error = function(ex) {})
  } else {
    setup(this, ..., verbose=less(verbose));
  }

  model <- getModel(this);

  # Generate layers?
  if (layers) {
    verbose && enter(verbose, "Generating image layers");

    # 1. Non-sample specific layers
    verbose && enter(verbose, "Non-sample specific layers");
    writeAxesLayers(this, chromosomes=chromosomes, zooms=zooms, ..., verbose=less(verbose));
    writeGridHorizontalLayers(this, chromosomes=chromosomes, zooms=zooms, ..., verbose=less(verbose));
    writeCytobandLayers(this, chromosomes=chromosomes, zooms=zooms, ..., verbose=less(verbose));
    verbose && exit(verbose);

    # 2. Sample specific layers
    verbose && enter(verbose, "Sample-specific layers");
    writeRawCopyNumberLayers(this, arrays=arrays, chromosomes=chromosomes, 
                                     zooms=zooms, ..., verbose=less(verbose));

    if (inherits(model, "CopyNumberSegmentationModel")) {
      writeCopyNumberRegionLayers(this, arrays=arrays, 
           chromosomes=chromosomes, zooms=zooms, ..., verbose=less(verbose));
    }
    verbose && exit(verbose);

    verbose && exit(verbose);
  } else {
    # Generate bitmap images
    writeGraphs(this, arrays=arrays, chromosomes=chromosomes, 
                                   zooms=zooms, ..., verbose=less(verbose));
  }

  # Update samples.js
  if (getParallelSafe(this)) {
    tryCatch({
      updateSamplesFile(this, ..., verbose=less(verbose));
    }, error = function(ex) {})
  } else {
    updateSamplesFile(this, ..., verbose=less(verbose));
  }

  # Write regions file
  if (inherits(model, "CopyNumberSegmentationModel")) {
    if (getParallelSafe(this)) {
      tryCatch({
        writeRegions(this, arrays=arrays, chromosomes=chromosomes, ..., verbose=less(verbose));
      }, error = function(ex) {})
    } else {
      writeRegions(this, arrays=arrays, chromosomes=chromosomes, ..., verbose=less(verbose));
    }
  }

  verbose && exit(verbose);
})


setMethodS3("display", "ChromosomeExplorer", function(this, filename="ChromosomeExplorer.html", ...) {
  NextMethod("display", this, filename=filename, ...);
})



##############################################################################
# HISTORY:
# 2009-05-19
# o Now testing for file permissions for creat-/writ-/updating files/dirs.
# 2008-12-13
# o Added argument 'zooms' to the constructor of ChromosomeExplorer.
#   Added methods get- and setZooms().
# 2008-06-05
# o Made updateSamplesFile(), writeAxesLayers(), writeGridHorizontalLayers(),
#   writeCytobandLayers(), writeRegions(), setup(), process() parallel safe.
# 2008-05-08
# o Now one can pass argument 'aliased=TRUE' to process() which causes the
#   ChromosomeExplorer and coupled CopyNumberSegmentationModel to return
#   tags that inferred from aliased full names.
# o Now updateSamplesFile() of ChromosomeExplorer passes arguments '...' 
#   to getFullNames().
# o Now getFullNames() of ChromosomeExplorer passes arguments '...' along.
# 2007-10-17
# o Now the  accepts CopyNumberChromosomalModel:s.ChromosomeExplorer
# o Added support for specifying the output version (at least for now).
# 2007-10-10
# o Added support for layers.  This should be backward compatible.
# 2007-10-09
# o Now process() writes cytoband/*.png layer images.
# o Added writeCytobandLayers() etc.
# 2007-09-30
# o Now the main HTML file for ChromosomeExplorer is ChromosomeExplorer.html,
#   which is in analogue to how the ArrayExplorer works.  This HTML now
#   loads to a different Javascript file with a different name so that already
#   existing index.html ChromosomeExplorer files will still work.
# 2007-09-04
# o Now ChromosomeExplorer recognizes CopyNumberSegmentationModel:s.
# 2007-05-08
# o Added default zoom levels to updateSamplesFile() for ChromosomeExplorer.  
#   This is applies the first time process() is called.
# 2007-03-19
# o Now ChromosomeExplorer extends Explorer.
# 2007-02-06
# o Now templates are in reports/templates/ and includes in reports/includes/.
# o Updated the path to <rootPath>/<dataSetName>/<tags>/<chipType>/<set>/.
# 2007-01-17
# o Now all 'arrays' arguments can contain array names.
# o Added getArrays() and setArrays() in order to focus on a subset of the
#   arrays in the model.
# 2007-01-16
# o Added getAlias() and setAlias(), and updated getName() accordingly.
#   This makes it easy to change the name of output set for subsets of
#   arrays.
# 2007-01-15
# o Added some more Rdoc comments.
# 2007-01-10
# o BUG FIX: setup() would only add index.html if includes/ were missing.
# 2007-01-08
# o Added display().
# 2007-01-07
# o TODO: The region filter for writeRegions() is hardwired for now.
# o Update to include regions.xls file on the CE output page too.
# o Now all HTML, CSS, and Javascript files are created too.
# 2007-01-04
# o Created.
##############################################################################
