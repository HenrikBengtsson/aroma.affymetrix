setMethodS3("stextChipType", "AffymetrixCdfFile", function(this, ...) {
  stextChipType(..., chipType=getChipType(this));
}, private=TRUE)



setMethodS3("getImage", "AffymetrixCdfFile", function(this, transforms=NULL, xrange=c(0,Inf), yrange=xrange, zrange=c(0,sqrt(2^16)), field=c("isPm"), levels=NULL, zoom=1, ..., verbose=FALSE) {
  require("EBImage") || throw("Package not loaded: EBImage.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'field':
#  field <- match.arg(field);
  
  # Argument 'zoom':
  zoom <- Arguments$getDouble(zoom, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Getting CDF image");

  verbose && enter(verbose, "Reading CDF image");
  fields <- unique(c("cell", field));
  data <- readDataFrame(this, fields=fields, ..., verbose=less(verbose,1));
  verbose && str(verbose, data);
  z <- vector(mode=mode(data[[field]]), 1);
  z <- matrix(z, nrow=nbrOfRows(this), ncol=nbrOfColumns(this));
  z[indexByRow(z, data[,"cell"])] <- data[[field]];
  rm(data);
  verbose && summary(verbose, as.vector(z));
  verbose && printf(verbose, "RAM: %.1fMB\n", object.size(z)/1024^2);
  verbose && exit(verbose);

  verbose && enter(verbose, "Transforming data");
  dim <- dim(z);
  mode <- mode(z);
  if (mode == "character") {
    if (is.null(levels)) {
      z <- factor(z);
    } else {
      z <- factor(z, levels=levels);
    }
    z <- as.integer(z);
  }
  dim(z) <- dim;
  verbose && str(verbose, z);
  verbose && exit(verbose);

  verbose && enter(verbose, "Creating Image");
  img <- getImage(z, scale=zoom, lim=zrange, ..., verbose=less(verbose, 1));
  verbose && print(verbose, img);
  verbose && exit(verbose);

  verbose && exit(verbose);

  # Return the 'field'
  attr(img, "field") <- field;

  img;
})


############################################################################
# HISTORY:
# 2006-09-16
# o Added getGenomeInformation() and stextChipType().
############################################################################
