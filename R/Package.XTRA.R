setMethodS3("getDefaultSettings", "Package", function(this, ...) {
  NULL;
})


setMethodS3("getSettings", "Package", function(this, ...) {
  # Use default settings
  settings <- getDefaultSettings(this);

  # Add package settings in options()
  optionName <- sprintf("%s.settings", getName(this));
  options <- getOption(optionName);

  if (length(options) > 0) {
    keys <- names(unlist(options));
    keys <- strsplit(keys, split=".", fixed=TRUE);

    for (key in keys) {
      settings[[key]] <- options[[key]];
    }
  }

  attr(settings, "option") <- optionName;

  settings;
})


setMethodS3("setSettings", "Package", function(this, value=list(), ...) {
  key <- sprintf("%s.settings", getName(this));
  args <- list(value);
  names(args) <- key;
  do.call("options", args);
})

setMethodS3("updateSettings", "Package", function(this, ...) {
  # Get all settings; if any is missing, the default value will be used
  settings <- getSettings(this);
  # Update the options list
  setSettings(this, settings);

  flat <- getLeaves(Options(settings));
  for (key in names(flat)) {
    # aromaSettings? /HB 2009-05-17  ('settings' does not work)
    setOption(aromaSettings, key, flat[[key]]);
  }
})

############################################################################
# HISTORY:
# 2009-02-22
# o Created.
############################################################################
