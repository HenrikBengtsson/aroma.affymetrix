setMethodS3("readGeneAssignments", "AffymetrixNetAffxCsvFile", function(this, ..., unique=TRUE, parse=TRUE, addNames=TRUE, flatten=TRUE, na.rm=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'unique':
  unique <- Arguments$getLogical(unique);

  # Argument 'unique':
  parse <- Arguments$getLogical(parse);

  # Argument 'addNames':
  addNames <- Arguments$getLogical(addNames);

  # Argument 'flatten':
  flatten <- Arguments$getLogical(flatten);

  # Argument 'na.rm':
  na.rm <- Arguments$getLogical(na.rm);
  na.rm <- (na.rm || flatten);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading gene assignments from ", class(this)[1]);
  verbose && print(verbose, this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  colClassPatterns <- c("*"="NULL", "probesetId|geneAssignment"="character");
  map <- readDataFrame(this, colClassPatterns=colClassPatterns, ...);
  verbose && cat(verbose, "Number of entries read: ", nrow(map));

  if (unique) {
    verbose && enter(verbose, "Dropping duplicated entries");

    n0 <- nrow(map);
    map <- unique(map);
    n <- nrow(map);
    verbose && printf(verbose, "Dropped %d (%.2f%%) out of %d duplicated entries\n", n-n0, 100*(n-n0)/n0, n0);
    verbose && cat(verbose, "Number of unique entries: ", nrow(map));

    verbose && exit(verbose);
  }

  res <- map$geneAssignment;
  names(res) <- map$probesetId;

  # Not used anymore
  rm(map);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Parse
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (parse) {
    verbose && enter(verbose, "Parsing entries");

    pairs <- strsplit(res, split=" /// ", fixed=TRUE);

    if (unique) {
      pairs <- lapply(pairs, FUN=unique);
    }

    ok <- !sapply(pairs, FUN=function(x) any(is.na(x)));
    if (na.rm) {
      pairs <- pairs[ok];
      ok <- !sapply(pairs, FUN=function(x) any(is.na(x)));
    }

    pairs[ok] <- lapply(pairs[ok], FUN=function(x) {
      x <- strsplit(x, split=" // ", fixed=TRUE);
      ns <- sapply(x, FUN=length);
      # Sanity check
      stopifnot(all(ns == 5));
      x <- unlist(x, use.names=FALSE);
      dimNA(x) <- c(5,NA);

      # In case there where duplicated subentries
      if (unique) {
        x <- unique(x);
      }

      t(x);
    });

    if (addNames) {
      # [1] HuGene-1_0-st-v1.na31.AFFX_README.NetAffx-CSV-Files.txt
      names <- c("accession", "geneSymbol", "geneTitle", "cytoband", "entrezGeneId");
      ok <- !sapply(pairs, FUN=function(x) any(is.na(x)));
      pairs[ok] <- lapply(pairs[ok], FUN=function(x) {
        colnames(x) <- names;
        x;
      });
    } # if (addNames)


    if (flatten) {
      verbose && enter(verbose, "Flattens list to table");
  
      verbose && enter(verbose, "Expanding unit names");
      ns <- sapply(pairs, FUN=NROW);
      unitNames <- mapply(names(ns), ns, FUN=rep);
      unitNames <- unlist(unitNames);
      # Sanity check
      stopifnot(length(unitNames) == sum(ns));
      verbose && exit(verbose);
  
      verbose && enter(verbose, "Stacking unit entries");
      data <- Reduce(rbind, pairs);
      rownames(data) <- NULL;
      verbose && exit(verbose);
  
      verbose && enter(verbose, "Building final table");
      # Sanity check
      stopifnot(length(unitNames) == nrow(data));
      data <- cbind(unitName=unitNames, data);
      rownames(data) <- NULL;
      verbose && exit(verbose);

      verbose && enter(verbose, "Coercing to data frame");
      data <- as.data.frame(data, stringsAsFactors=FALSE);
      verbose && exit(verbose);
  
      pairs <- data;
  
      verbose && exit(verbose);
    } # if (flatten)

    res <- pairs;

    verbose && exit(verbose);
  } # if (parse)

  verbose && exit(verbose);

  res;
})  # readGeneAssignments()



##############################################################################
# HISTORY:
# 2011-04-06
# o Added readGeneAssignments() for AffymetrixNetAffxCsvFile.
# o Created.
##############################################################################
