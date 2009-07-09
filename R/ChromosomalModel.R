###########################################################################/**
# @RdocClass ChromosomalModel
#
# @title "The ChromosomalModel class"
#
# \description{
#  @classhierarchy
#
#  This \emph{abstract} class represents a chromosomal model.
# }
# 
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "ChipEffectSetTuple".}
#   \item{tags}{A @character @vector of tags.}
#   \item{genome}{A @character string specifying what genome is process.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   This class requires genome information annotation files for 
#   every chip type.
# }
#
# @author
#*/###########################################################################
setConstructorS3("ChromosomalModel", function(cesTuple=NULL, tags="*", genome="Human", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cesTuple':
  if (!is.null(cesTuple)) {
    # Coerce to ChipEffectSetTuple
    if (!inherits(cesTuple, "ChipEffectSetTuple")) {
      cesTuple <- ChipEffectSetTuple(cesTuple);
    }

    # Assert that we are dealing with CnChipEffectSet:s
    for (ces in getListOfSets(cesTuple)) {
      # Assert special properties for CnChipEffectSet:s AD HOC /HB 2006-12-20
      if (inherits(ces, "CnChipEffectSet")) {
        # Currently only total copy-number estimates are accepted
        if (!ces$combineAlleles) {
          throw("Unsupported copy-number chip effects. Currently only total copy-number estimates are supported: ces$combineAlleles == FALSE");
        }
      }
    }
  }

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }
 

  this <- extend(Object(), "ChromosomalModel",
    .alias = NULL,
    .cesTuple = cesTuple,
    .chromosomes = NULL,
    .tags = tags,
    .genome = genome
  );

  # Validate?
  if (!is.null(this$.cesTuple)) {
    # Validate genome
    pathname <- getGenomeFile(this);
  }

  this;
}, abstract=TRUE)


setMethodS3("as.character", "ChromosomalModel", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", getTags(this, collapse=",")));
  s <- c(s, paste("Chip type (virtual):", getChipType(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  nbrOfChipTypes <- nbrOfChipTypes(this);
  s <- c(s, sprintf("Number of chip types: %d", nbrOfChipTypes));
  s <- c(s, "Chip-effect set:");
  cesList <- getListOfSets(getSetTuple(this));
  chipTypes <- sapply(cesList, FUN=function(ces) {
    getChipType(ces);
  })
  for (kk in seq(along=cesList)) {
    s <- c(s, sprintf("Chip type #%d of %d ('%s'):", 
                                    kk, nbrOfChipTypes, chipTypes[kk]));
    s <- c(s, "Chip-effect set:");
    ces <- cesList[[kk]];
    s <- c(s, as.character(ces));
  }
#  s <- c(s, "Genome information:", as.character(getGenomeInformation(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, protected=TRUE)


setMethodS3("clearCache", "ChromosomalModel", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c()) {
    this[[ff]] <- NULL;
  }

  if (!is.null(this$.cesTuple))
    clearCache(this$.cesTuple);

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
})


setMethodS3("getRootPath", "ChromosomalModel", function(this, ...) {
  tag <- getAsteriskTags(this)[1];
  sprintf("%sData", tolower(tag));
})


setMethodS3("getParentPath", "ChromosomalModel", function(this, ...) {
  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # The full path
  path <- filePath(rootPath, fullname, expandLinks="any");

  # Create path?
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create directory: ", path);
  }

  path;
})

setMethodS3("getPath", "ChromosomalModel", function(this, ...) {
  path <- getParentPath(this, ...);

  # Chip type
  chipType <- getChipType(this);

	  # The full path
  path <- filePath(path, chipType, expandLinks="any");

  # Create path?
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})

setMethodS3("getReportPath", "ChromosomalModel", function(this, ...) {
  rootPath <- "reports";

  # Data set name
  name <- getName(this);

  # Data set tags
  tags <- getTags(this, collapse=",");

  # Get chip type
  chipType <- getChipType(this);

  # Image set
  set <- getSetTag(this);

  # The report path
  path <- filePath(rootPath, name, tags, chipType, set, expandLinks="any");

  path;
}, protected=TRUE)



setMethodS3("getSetTuple", "ChromosomalModel", function(this, ...) {
  this$.cesTuple;
}, protected=TRUE)


setMethodS3("getListOfChipEffectSets", "ChromosomalModel", function(this, ...) {
  cesTuple <- getSetTuple(this);
  getListOfSets(cesTuple);
})



###########################################################################/**
# @RdocMethod nbrOfChipTypes
#
# @title "Gets the number of chip types"
#
# \description{
#  @get "title" used in the model.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfChipTypes", "ChromosomalModel", function(this, ...) {
  nbrOfChipTypes(getSetTuple(this), ...);
})



setMethodS3("getListOfUnitNamesFiles", "ChromosomalModel", function(this, ...) {
  getListOfUnitNamesFiles(getSetTuple(this), ...);
}, private=TRUE)

setMethodS3("getListOfUnitTypesFiles", "ChromosomalModel", function(this, ...) {
  getListOfUnitTypesFiles(getSetTuple(this), ...);
}, private=TRUE)



setMethodS3("getChipTypes", "ChromosomalModel", function(this, ...) {
  getChipTypes(getSetTuple(this), ...);
})


###########################################################################/**
# @RdocMethod getChipType
#
# @title "Gets a label for all chip types merged"
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
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipType", "ChromosomalModel", function(this, ...) {
  getChipTypes(this, merge=TRUE, ...);
})


###########################################################################/**
# @RdocMethod getTableOfArrays
#
# @title "Gets a table of arrays"
#
# \description{
#  @get "title" showing their availability across chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a \eqn{NxK} @matrix of @integers where \eqn{N} is the total number 
#  of arrays and \eqn{K} is the number of chip types in the model.  The row 
#  names are the names of the arrays, and the column names are the chip types.
#  If data is available for array \eqn{n} and chip type \eqn{k}, cell 
#  \eqn{(n,k)} has value \eqn{n}, otherwise @NA.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTableOfArrays", "ChromosomalModel", function(this, ...) {
  getTableOfArrays(getSetTuple(this), ...);
})


setMethodS3("getNames", "ChromosomalModel", function(this, ...) {
  rownames(getTableOfArrays(this, ...));
})


setMethodS3("getFullNames", "ChromosomalModel", function(this, ...) {
  getFullNames(getSetTuple(this), ...);
})



###########################################################################/**
# @RdocMethod getArrays
#
# @title "Gets the names of the arrays"
#
# \description{
#  @get "title" available in the model.
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
setMethodS3("getArrays", "ChromosomalModel", function(this, ...) {
  getNames(this, ...);
})




###########################################################################/**
# @RdocMethod indexOfArrays
#
# @title "Gets the indices of the arrays"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @character @vector of arrays names.
#     If @NULL, all arrays are considered.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("indexOfArrays", "ChromosomalModel", function(this, arrays=NULL, ...) {
  indexOfArrays(getSetTuple(this), arrays=arrays, ...);
}, private=TRUE)



###########################################################################/**
# @RdocMethod nbrOfArrays
#
# @title "Gets the number of arrays"
#
# \description{
#  @get "title" used in the model.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfArrays", "ChromosomalModel", function(this, ...) {
  length(getNames(this, ...));
})


setMethodS3("getName", "ChromosomalModel", function(this, collapse="+", ...) {
  name <- getAlias(this);

  if (is.null(name))
    name <- getName(getSetTuple(this), ...);

  name;
})

setMethodS3("getAlias", "ChromosomalModel", function(this, ...) {
  this$.alias;
})


setMethodS3("setAlias", "ChromosomalModel", function(this, alias=NULL, ...) {
  # Argument 'alias':
  if (!is.null(alias)) {
    alias <- Arguments$getCharacter(alias, length=c(1,1));
  }
  this$.alias <- alias;
  invisible(this);
})





setMethodS3("getAsteriskTags", "ChromosomalModel", function(this, collapse=NULL, ...) {
  # Create a default asterisk tags for any class by extracting all
  # capital letters and pasting them together, e.g. AbcDefGhi => ADG.
  name <- class(this)[1];

  # Remove any 'Model' suffixes
  name <- gsub("Model$", "", name);

  name <- capitalize(name);

  # Vectorize
  name <- strsplit(name, split="")[[1]];

  # Identify upper case
  name <- name[(toupper(name) == name)];

  # Paste
  name <- paste(name, collapse="");

  tag <- name;
}, protected=TRUE)




setMethodS3("getTags", "ChromosomalModel", function(this, collapse=NULL, ...) {
  tags <- getTags(getSetTuple(this), collapse=collapse, ...);

  # Add model tags
  tags <- c(tags, this$.tags);

  # Update default tags
  asteriskTags <- getAsteriskTags(this, collapse=",");
  if (length(asteriskTags) == 0)
    asteriskTags <- "";
  tags[tags == "*"] <- asteriskTags;
  tags <- tags[nchar(tags) > 0];
  tags <- unlist(strsplit(tags, split=","), use.names=TRUE);

  # Keep non-empty tags
  tags <- tags[nchar(tags) > 0];

  # Get unique tags
  tags <- locallyUnique(tags);

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("getFullName", "ChromosomalModel", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
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
setMethodS3("getChromosomes", "ChromosomalModel", function(this, ...) {
  gis <- getListOfGenomeInformations(this);
  chromosomes <- lapply(gis, getChromosomes);
  chromosomes <- unlist(chromosomes, use.names=TRUE);
  chromosomes <- sort(unique(chromosomes));
  chromosomes;
})




setMethodS3("getListOfAromaUgpFiles", "ChromosomalModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving UGP files");
  unfList <- getListOfUnitNamesFiles(this);
  ugpList <- lapply(unfList, getAromaUgpFile, verbose=less(verbose));
  verbose && exit(verbose);

  ugpList;
})


setMethodS3("getChipEffectFiles", "ChromosomalModel", function(this, ...) {
  setTuple <- getSetTuple(this);
  getArrayTuple(setTuple, ...);
})




setMethodS3("getGenome", "ChromosomalModel", function(this, ...) {
  this$.genome;
})


setMethodS3("getGenomeFile", "ChromosomalModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  fullname <- getGenome(this);
  pattern <- sprintf("^%s.*,chromosomes.txt$", fullname);

  # 1. Search in the regular places
  pathname <- findAnnotationData(name=fullname, set="genomes", 
                            pattern=pattern, ..., verbose=less(verbose, 10));

  # 2. As a backup, search in the <pkg>/annotationData/ directory
  if (is.null(pathname)) {
    verbose && enter(verbose, "Search among package's annotationData/");
    path <- system.file("annotationData", package="aroma.affymetrix");
    verbose && cat(verbose, "Path: ", path);
    pathname <- findAnnotationData(name=fullname, set="genomes", 
                pattern=pattern, ..., paths=path, verbose=less(verbose, 10));
    verbose && exit(verbose);
  }

  if (is.null(pathname)) {
    throw("Failed to locate a genome annotation data file: ", fullname);
  }

  pathname;
}, protected=TRUE)


setMethodS3("setGenome", "ChromosomalModel", function(this, genome, tags=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genome':
  genome <- Arguments$getCharacter(genome, length=c(1,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  oldGenome <- this$.genome;

  fullname <- paste(c(genome, tags), collapse=",");
  verbose && cat(verbose, "Fullname: ", fullname);

  # Verify that there is an existing genome file
  tryCatch({
    this$.genome <- fullname;
    pathname <- getGenomeFile(this, verbose=less(verbose, 10));
  }, error = function(ex) {
    this$.genome <- oldGenome;
    throw(ex$message);
  })

  invisible(oldGenome);
})



setMethodS3("getGenomeData", "ChromosomalModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading genome chromosome annotation file");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get genome annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Searching for the file");
  # Search annotationData/genomes/
  pathname <- getGenomeFile(this, verbose=less(verbose, 10));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading data file");
  verbose && cat(verbose, "Pathname: ", pathname);
  data <- readTable(pathname, header=TRUE, 
                            colClasses=c(nbrOfBases="integer"), row.names=1);
  verbose && exit(verbose);

  verbose && enter(verbose, "Translating chromosome names");
  chromosomes <- row.names(data);
  map <- c("X"=23, "Y"=24, "Z"=25);
  for (kk in seq(along=map)) {
    chromosomes <- gsub(names(map)[kk], map[kk], chromosomes, fixed=TRUE);
  }
  row.names(data) <- chromosomes;
  verbose && exit(verbose);

  verbose && exit(verbose);

  data;
}, protected=TRUE)


setMethodS3("fit", "ChromosomalModel", abstract=TRUE);


setMethodS3("getSetTag", "ChromosomalModel", function(this, ...) {
  tolower(getAsteriskTags(this)[1]);
}, private=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# BEGIN: AFFX specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setMethodS3("getListOfGenomeInformations", "ChromosomalModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving genome informations");
  cdfList <- getListOfCdfs(getSetTuple(this), ...);
  giList <- lapply(cdfList, getGenomeInformation, verbose=less(verbose));
  verbose && exit(verbose);

  giList;
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# END: AFFX specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



##############################################################################
# HISTORY:
# 2009-07-08
# o Added getListOfUnitTypesFiles() for ChromosomalModel.
# 2009-01-26
# o Removed get[]ListOfCdfs() from ChromosomalModel.
# o Removed deprectated get[]ListOfChipEffects() from ChromosomalModel.
# o Added getListOfAromaUgpFiles() to ChromosomalModel.
# o Added getListOfUnitNamesFiles() to ChromosomalModel.
# 2007-09-25
# o Extracted ChromosomalModel from CopyNumberSegmentationModel.  For 
#   previous HISTORY, see that class.
##############################################################################
