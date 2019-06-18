# Repitools

<details>

* Version: 1.30.0
* Source code: https://github.com/cran/Repitools
* Date/Publication: 2019-05-02
* Number of recursive dependencies: 118

Run `revdep_details(,"Repitools")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    ...
    writeWig,AffymetrixCelSet: no visible global function definition for
      ‘extractMatrix’
    Undefined global functions or variables:
      Arguments AromaCellCpgFile AromaCellPositionFile
      AromaCellSequenceFile DNAString abline axis barplot bxp countBases
      dbeta dev.off embed enter exit extract extractMatrix filter getCdf
      getCellIndices getChipType getMainCdf getNames grid indexOf kmeans
      layout legend lines lm lowess matchPattern matlines matplot mtext
      nbrOfArrays nbrOfUnits p.adjust par pdf persp plot plot.new
      plot.window points polygon popState predict pt pushState qnorm
      rainbow read.table rect str t.test text title verbose
    Consider adding
      importFrom("grDevices", "dev.off", "pdf", "rainbow")
      importFrom("graphics", "abline", "axis", "barplot", "bxp", "grid",
                 "layout", "legend", "lines", "matlines", "matplot", "mtext",
                 "par", "persp", "plot", "plot.new", "plot.window", "points",
                 "polygon", "rect", "text", "title")
      importFrom("stats", "dbeta", "embed", "filter", "kmeans", "lm",
                 "lowess", "p.adjust", "predict", "pt", "qnorm", "t.test")
      importFrom("utils", "read.table", "str")
    to your NAMESPACE file.
    ```

# TIN

<details>

* Version: 1.16.0
* Source code: https://github.com/cran/TIN
* Date/Publication: 2019-05-02
* Number of recursive dependencies: 114

Run `revdep_details(,"TIN")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    ...
    scatterPlot: no visible global function definition for ‘pdf’
    scatterPlot: no visible global function definition for ‘bmp’
    scatterPlot: no visible global function definition for ‘plot’
    scatterPlot: no visible global function definition for ‘ave’
    scatterPlot: no visible global function definition for ‘axis’
    scatterPlot: no visible global function definition for ‘text’
    scatterPlot: no visible global function definition for ‘mtext’
    scatterPlot: no visible global function definition for ‘points’
    scatterPlot: no visible global function definition for ‘dev.off’
    Undefined global functions or variables:
      ave axis bmp colorRampPalette data dev.off dist hclust hist jpeg
      median mtext par pdf plot png points postscript quantile read.table
      text
    Consider adding
      importFrom("grDevices", "bmp", "colorRampPalette", "dev.off", "jpeg",
                 "pdf", "png", "postscript")
      importFrom("graphics", "axis", "hist", "mtext", "par", "plot",
                 "points", "text")
      importFrom("stats", "ave", "dist", "hclust", "median", "quantile")
      importFrom("utils", "data", "read.table")
    to your NAMESPACE file.
    ```

