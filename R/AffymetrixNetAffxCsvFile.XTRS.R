setMethodS3("readGeneAssignments", "AffymetrixNetAffxCsvFile", function(this, ..., unique=TRUE, parse=TRUE, fields=NULL, flatten=TRUE, na.rm=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'unique':
  unique <- Arguments$getLogical(unique);

  # Argument 'unique':
  parse <- Arguments$getLogical(parse);

  knownFields <- c("accession", "geneSymbol", "geneTitle", "cytoband", "entrezGeneId");
  if (!is.null(fields)) {
    fields <- match.arg(fields, choices=knownFields, several.ok=TRUE);
  }

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

    
    nFields <- c(2L,5L);

    # Infer number of fields from data?
    if (length(nFields) > 1) {
      x <- pairs[ok][[1]];
      x <- strsplit(x, split=" // ", fixed=TRUE);
      ns <- sapply(x, FUN=length);
      n <- unique(ns);
      # Sanity check
      stopifnot(length(n) == 1);
      stopifnot(is.element(n, nFields));
      nFields <- n;
    }
    knownFields <- knownFields[1:nFields];

    pairs[ok] <- lapply(pairs[ok], FUN=function(x) {
      x <- strsplit(x, split=" // ", fixed=TRUE);
      ns <- sapply(x, FUN=length);
      # Sanity check
      stopifnot(all(ns == nFields));
      x <- unlist(x, use.names=FALSE);
      dimNA(x) <- c(nFields,NA);

      x <- t(x);

      # In case there where duplicated subentries
      if (unique) {
        x <- unique(x);
      }

      x;
    });


    verbose && enter(verbose, "Adding column names");
    verbose && cat(verbose, "Column names: ", hpaste(knownFields));

    # [1] HuGene-1_0-st-v1.na31.AFFX_README.NetAffx-CSV-Files.txt
    ok <- !sapply(pairs, FUN=function(x) any(is.na(x)));
    pairs[ok] <- lapply(pairs[ok], FUN=function(x) {
      colnames(x) <- knownFields;
      x;
    });

    verbose && exit(verbose);


    # Check argument 'fields' again
    if (is.null(fields)) {
      fields <- knownFields;
    }

    verbose && cat(verbose, "Fields: ", hpaste(fields));

    ok <- !sapply(pairs, FUN=function(x) any(is.na(x)));
    pairs[ok] <- lapply(pairs[ok], FUN=function(x) {
      x[,fields,drop=FALSE];
    });

    if (unique) {
      pairs[ok] <- lapply(pairs[ok], FUN=unique);
    }

    if (flatten) {
      verbose && enter(verbose, "Flattens list to table");

      verbose && enter(verbose, "Identifying blocks of unique sizes");
      ns <- sapply(pairs, FUN=NROW);
      uns <- unique(ns);
      verbose && print(verbose, uns);
      verbose && exit(verbose);

      verbose && enter(verbose, "Stacking by block size");
      unitNames <- idxs <- c();
      for (ii in seq(along=uns)) {
        verbose && enter(verbose, sprintf("Size %d of %d", ii, length(uns)));
        n <- uns[ii];
        keep <- which(ns == n);
        verbose && str(verbose, keep);

        pairsII <- pairs[keep];
        unitNamesII <- rep(names(pairsII), each=n);
        idxsII <- rep(seq.int(n), times=length(pairsII));

        unitNames <- c(unitNames, unitNamesII);
        idxs <- c(idxs, idxsII);

        verbose && exit(verbose);
      } # for (ii ...)
  
      # Sanity checks
      stopifnot(length(unitNames) == sum(ns));
      stopifnot(length(idxs) == sum(ns));

      verbose && exit(verbose);
  
      verbose && enter(verbose, "Stacking unit entries");
      data <- Reduce(rbind, pairs);
      rownames(data) <- NULL;
      verbose && exit(verbose);
  
      verbose && enter(verbose, "Building final table");
      # Sanity check
      stopifnot(length(unitNames) == nrow(data));
      data <- cbind(unitName=unitNames, index=idxs, data);
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
# 2011-04-07
# o Now readGeneAssignments() for AffymetrixNetAffxCsvFile handles both
#   *.probeset.csv (2 fields) and *.transcript.csv (5 fields), at least
#   for the HuGene-1_0-st-v1 chip type.
# 2011-04-06
# o Added readGeneAssignments() for AffymetrixNetAffxCsvFile.
# o Created.
##############################################################################
