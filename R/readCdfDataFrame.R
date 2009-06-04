#########################################################################/**
# @RdocFunction readCdfDataFrame
#
# @title "Reads units (probesets) from an Affymetrix CDF file"
#
# @synopsis 
# 
# \description{
#  @get "title". Gets all or a subset of units (probesets).
# }
# 
# \arguments{
#  \item{filename}{The filename of the CDF file.}
#  \item{units}{An @integer @vector of unit indices
#    specifying which units to be read.  If @NULL, all are read.}
#  \item{groups}{An @integer @vector of group indices
#    specifying which groups to be read.  If @NULL, all are read.}
#  \item{cells}{An @integer @vector of cell indices
#    specifying which cells to be read.  If @NULL, all are read.}
#  \item{fields}{A @character @vector specifying what fields to read.
#    If @NULL, all unit, group and cell fields are returned.}
#  \item{drop}{If @TRUE and only one field is read, then a @vector
#    (rather than a single-column @data.frame) is returned.}
#  \item{verbose}{An @integer specifying the verbose level. If 0, the
#    file is parsed quietly.  The higher numbers, the more details.}
# }
# 
# \value{
#   An NxK @data.frame or a @vector of length N.
# }
#
# \author{Henrik Bengtsson}
# 
# @examples "../incl/readCdfDataFrame.Rex"
# 
# \seealso{
#   For retrieving the CDF as a @list structure, see 
#   @see "affxparser::readCdfUnits".
# }
# 
# \references{
#   [1] Affymetrix Inc, Affymetrix GCOS 1.x compatible file formats,
#       June 14, 2005.
#       \url{http://www.affymetrix.com/support/developer/}
# }
#
# @keyword "file"
# @keyword "IO"
# @keyword "internal"
#*/######################################################################### 
readCdfDataFrame <- function(filename, units=NULL, groups=NULL, cells=NULL, fields=NULL, drop=TRUE, verbose=0) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'filename':
  filename <- file.path(dirname(filename), basename(filename));
  if (!file.exists(filename))
    stop("File not found: ", filename);

  # Argument 'units':
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- as.integer(units);
    if (any(units < 1))
      stop("Argument 'units' contains non-positive indices.");
  } else {
    stop("Argument 'units' must be numeric or NULL: ", class(units)[1]);
  }

  # Argument 'groups':
  if (is.null(groups)) {
  } else if (is.numeric(groups)) {
    groups <- as.integer(groups);
    if (any(groups < 1))
      stop("Argument 'groups' contains non-positive indices.");
  } else {
    stop("Argument 'groups' must be numeric or NULL: ", class(groups)[1]);
  }

  # Argument 'fields':
##   knownUnitFields <- c("unit", "unitName", "unitDirection", "nbrOfUnitAtoms", "unitSize", "unitNumber", "unitType", "nbrOfGroups", "mutationType");
##   knownGroupFields <- c("group", "groupName", "nbrOfGroupAtoms", "groupSize", "firstAtom", "lastAtom", "groupDirection");
##   knownCellFields <- c("cell", "x", "y", "probeSequence", "feat", "qual", "expos", "pos", "cbase", "pbase", "tbase", "atom", "index");

  if (is.null(fields)) {
    knownUnitFields <- c("unit", "unitName", "unitType", 
                         "unitDirection", "unitNbrOfAtoms");
    knownGroupFields <- c("group", "groupName", "groupDirection",
                          "groupNbrOfAtoms");
    knownCellFields <- c("cell", "x", "y", "pbase", "tbase", 
                         "indexPos", "atom", "expos");
    fields <- c(knownUnitFields, knownGroupFields, knownCellFields);
  }
  
  # Argument 'verbose':
  if (length(verbose) != 1)
    stop("Argument 'verbose' must be a single integer.");
  verbose <- as.integer(verbose);
  if (!is.finite(verbose))
    stop("Argument 'verbose' must be an integer: ", verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Prepare the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  readFields <- c(fields, "cell");  # Need to read one cell field!
  readFields <- unique(readFields);
  
  # Unit fields
  readUnitType <- ("unitType" %in% readFields);
  readUnitDirection <- ("unitDirection" %in% readFields);
  readUnitNumber <- ("unitNumber" %in% readFields);
  readUnitAtomNumbers <- ("unitNbrOfAtoms" %in% readFields);

  # Group fields
  readGroupDirection <- ("groupDirection" %in% readFields);
  readGroupAtomNumbers <- ("groupNbrOfAtoms" %in% readFields);

  # Cell fields
  readXY <- any(c("x", "y") %in% readFields);
  readIndices <- ("cell" %in% readFields);
  readBases <- any(c("tbase", "pbase") %in% readFields);
  readIndexpos <- ("indexPos" %in% readFields);
  readExpos <- ("expos" %in% readFields);
  readAtoms <- ("atom" %in% readFields);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Query the CDF
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdf <- readCdf(filename, units=units, readXY=readXY, readBases=readBases,
    readIndexpos=readIndexpos, readAtoms=readAtoms, 
    readUnitType=readUnitType, readUnitDirection=readUnitDirection, 
    readUnitNumber=readUnitNumber, readUnitAtomNumbers=readUnitAtomNumbers,
    readGroupAtomNumbers=readGroupAtomNumbers, 
    readGroupDirection=readGroupDirection, readIndices=readIndices, 
    verbose=verbose-1);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Flatten CDF list structure unit by unit
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(units))
    units <- seq(along=cdf);
  groupIdxs <- groups;

  unitNames <- names(cdf);
  for (uu in seq(along=cdf)) {
    unit <- .subset2(cdf, uu);
    unitName <- .subset(unitNames, uu);

    if (verbose >= 1) {
      if (uu %% 500 == 1) {
        unitsLeft <- length(cdf) - uu + 1;
        cat(unitsLeft, ", ", sep="");
      }
    }

    groups <- .subset2(unit, "groups");
    unit[["groups"]] <- NULL;

    # Translate unit names (has to be done here because not unique)
    names <- names(unit);
    names <- sub("ncells", "unitNbrOfCells", names); 
    names <- sub("natoms", "unitNbrOfAtoms", names); 
    names <- sub("unitnumber", "unitNumber", names); 
    names(unit) <- names;

    unitData <- list(unit=.subset(units, uu), unitName=unitName);
    unitData <- c(unitData, unit);

    # Extract groups of interest?
    if (is.null(groupIdxs)) {
      ggs <- seq(along=groups);
    } else {
      keep <- which(seq(along=groups) %in% groupIdxs);
      groups <- .subset(groups, keep);
      ggs <- groupIdxs;
    }

    # Flatten (group, cell) data
    groupNames <- names(groups);
    for (gg in seq(along=ggs)) {
      group <- .subset2(groups, gg);
      groupName <- .subset(groupNames, gg);

      groupData <- list(group=.subset(ggs, gg), groupName=groupName);

      # Extract group fields
      keys <- c("groupdirection", "natoms", "ncellsperatom", "ncells");
      idxs <- which(names(group) %in% keys);
      if (length(idxs) > 0) {
        groupData <- c(groupData, .subset(group, idxs));
        group <- .subset(group, -idxs);
      }

      # Extract cell fields
      cellData <- as.data.frame(group, stringsAsFactors=FALSE);

      # Extract cells of interest?
      if (!is.null(cells)) {
        keep <- (seq(length=nrow(cellData)) %in% cells);
        cellData <- cellData[keep,,drop=FALSE];
#        rm(keep);
      }

      # Expand group fields
      nbrOfCells <- nrow(cellData);
      for (key in names(groupData)) {
        groupData[[key]] <- rep(.subset2(groupData, key), times=nbrOfCells);
      }
      groupData <- as.data.frame(groupData, stringsAsFactors=FALSE);

      group <- cbind(groupData, cellData);
#      rm(groupData, cellData);

      groups[[gg]] <- group;
    }

    # Stack (rbind) groups
    stackedGroups <- NULL;
    for (gg in seq(along=groups)) {
      stackedGroups <- rbind(stackedGroups, .subset2(groups, gg));
    }
#    rm(groups);

    nbrOfCells <- nrow(stackedGroups);
    for (key in names(unitData)) {
      unitData[[key]] <- rep(.subset2(unitData, key), times=nbrOfCells);
    }
    unitData <- as.data.frame(unitData, stringsAsFactors=FALSE);

    unit <- cbind(unitData, stackedGroups);
#    rm(unitData, stackedGroups);

    cdf[[uu]] <- unit;
  } # for (uu ...)

  if (verbose >= 1) {
    cat("0.\n");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Flatten the remaining list structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Allocate "data frame" of right size
  unitSizes <- sapply(cdf, FUN=nrow);
  nbrOfCells <- sum(unitSizes);
  dataTypes <- sapply(.subset2(cdf, 1), FUN=storage.mode);
  df <- vector("list", length(dataTypes));
  names(df) <- names(dataTypes);
  for (key in names(df)) {
    df[[key]] <- vector(.subset2(dataTypes, key), length=nbrOfCells);
  }

  # Copy values from the CDF list structure
  offset <- 0;
  for (uu in seq(along=cdf)) {
    data <- .subset2(cdf, uu);
    nrow <- nrow(data);
    idxs <- offset + 1:nrow;
    for (key in names(df)) {
      df[[key]][idxs] <- .subset2(data, key);
    }
    offset <- offset + nrow;
    cdf[[uu]] <- NA;
  }

  names <- names(df);

  # Translate unit names
  names <- sub("unittype", "unitType", names); 
  names <- sub("unitdirection", "unitDirection", names); 
  names <- sub("ncellsperatom", "unitNbrOfCellsPerAtom", names); 
  # Translate group names
  names <- sub("groupdirection", "groupDirection", names); 
  names <- sub("natoms", "groupNbrOfAtoms", names); 
  names <- sub("ncellsperatom", "groupNbrOfCellsPerAtom", names); 
  names <- sub("ncells", "groupNbrOfCells", names); 
  # Translate cell names
  names <- sub("indices", "cell", names); 
  names <- sub("indexpos", "indexPos", names); 
  names(df) <- names;

  # Extract fields of interest
  unknown <- setdiff(fields, names(df))
  if (length(unknown) > 0) {
    warning("Some of the fields were not read: ", 
                                   paste(unknown, collapse=", "));
  }
  fields <- intersect(fields, names(df));
  df <- .subset(df, fields);


  # Make it a valid data frame
  if (drop && length(df) == 1) {
    df <- .subset2(df, 1);
  } else {
    attr(df, "row.names") <- .set_row_names(nbrOfCells);
    attr(df, "class") <- "data.frame";
  }

  df;
} # readCdfDataFrame()





# TO DO:
.readCdfDataFrameFinal <- function(filename, units=NULL, groups=NULL, unitFields=NULL, groupFields=NULL, cellFields=NULL, verbose=0) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'filename':
  filename <- file.path(dirname(filename), basename(filename));
  if (!file.exists(filename))
    stop("File not found: ", filename);

  # Argument 'units':
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- as.integer(units);
    if (any(units < 1))
      stop("Argument 'units' contains non-positive indices.");
  } else {
    stop("Argument 'units' must be numeric or NULL: ", class(units)[1]);
  }

  # Argument 'groups':
  if (is.null(groups)) {
  } else if (is.numeric(groups)) {
    groups <- as.integer(groups);
    if (any(groups < 1))
      stop("Argument 'groups' contains non-positive indices.");
  } else {
    stop("Argument 'groups' must be numeric or NULL: ", class(groups)[1]);
  }

  # Argument 'unitFields':
  knownFields <- c("unit", "unitName", "unitDirection", "nbrOfUnitAtoms", "unitSize", "unitNumber", "unitType", "nbrOfGroups", "mutationType");
  
  # Argument 'groupFields':
  knownFields <- c("group", "groupName", "nbrOfGroupAtoms", "groupSize", "firstAtom", "lastAtom", "groupDirection");

  # Argument 'cellFields':
  knownFields <- c("cell", "x", "y", "probeSequence", "feat", "qual", "expos", "pos", "cbase", "pbase", "tbase", "atom", "index");


  # Argument 'verbose':
  if (length(verbose) != 1)
    stop("Argument 'verbose' must be a single integer.");
  verbose <- as.integer(verbose);
  if (!is.finite(verbose))
    stop("Argument 'verbose' must be an integer: ", verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Call the C function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
##  df <- .Call("R_affx_get_cdf_dataframe", filename, units, 
##               unitFields, groupFields, cellFields,
##               unitNames, groupNames,
##               verbose, PACKAGE="affxparser");

  df;
} # .readCdfDataFrameFinal()



############################################################################
# HISTORY:
# 2008-04-03 [HB]
# o Now the renaming of fields is mostly done at the end and not in every
#   iteration.  That speeds up the process 5-10%.
# o Replaced all [[() and [() with .subset2() and .subset(), respectively.
#   This should make it a little bit faster.
# o For my record, some benchmarking on readCdfDataFrame():
#   Mapping10K_Xba142.cdf: 408,400x16, 23.3Mb, 610 seconds.
# 2008-03-24 [HB]
# o Created a readCdfDataFrame() that uses readCdf() to read the data and
#   then restructures it.  This will be used as a reference for the final
#   implementation.
# o Created a stub for readCdfDataFrame().
############################################################################  
