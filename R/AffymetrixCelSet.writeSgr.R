
setMethodS3("writeSgr", "AffymetrixCelSet", function(this, units=NULL, ..., tags=getTags(this), verbose=FALSE, fileExtension="sgr", fileSep="\t",nbrOfFigures=3) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  
  chrText <- paste(c(1:22,"X","Y","M"),sep="")
  names(chrText) <- 1:25
  
  nbrOfArrays <- nbrOfArrays(this)
  cdf <- getCdf(this)

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Reading indices");
  indices <- getCellIndices(cdf,units=units,stratifyBy="pm")
  indices <- unlist(indices,use.names=FALSE)
  verbose && exit(verbose);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read probe genomic location
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Gathering probe genomic locations");
  acp <- AromaCellPositionFile$byChipType(getChipType(cdf))
  
  ch <- paste("chr", chrText[ as.character(acp[indices,1,drop=TRUE]) ], sep="")
  pos <- acp[indices,2,drop=TRUE]
  verbose && exit(verbose);
  
  sampleNames <- getNames(this)
  
  for(ii in seq_len(nbrOfArrays)) {
  
    verbose && enter(verbose, sprintf("Gathering and writing data for ", sampleNames[ii]));
    cf <- getFile(this, ii)
	
	data <- extractMatrix(cf, cells=indices, verbose=verbose)
	data <- log2(data)
	
	fileName <- paste(sampleNames[ii],paste(tags,collapse=","),fileExtension,sep=".")
	
	write.table( cbind(ch,pos,round(data,nbrOfFigures)), fileName, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=fileSep)
    verbose && exit(verbose);
  
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
})

