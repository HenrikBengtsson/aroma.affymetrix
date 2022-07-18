# TIN

<details>

* Version: 1.28.0
* GitHub: NA
* Source code: https://github.com/cran/TIN
* Date/Publication: 2022-04-26
* Number of recursive dependencies: 124

Run `revdep_details(, "TIN")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    aberrantExonUsage: no visible global function definition for ‘quantile’
    aberrantExonUsage: no visible global function definition for ‘ave’
    clusterPlot: no visible global function definition for ‘dist’
    clusterPlot: no visible global function definition for ‘hclust’
    clusterPlot: no visible global function definition for
      ‘colorRampPalette’
    clusterPlot: no visible global function definition for ‘par’
    clusterPlot: no visible global function definition for ‘png’
    clusterPlot: no visible global function definition for ‘jpeg’
    clusterPlot: no visible global function definition for ‘postscript’
    ...
      ave axis bmp colorRampPalette data dev.off dist hclust hist jpeg
      median mtext par pdf png points postscript quantile read.table text
    Consider adding
      importFrom("grDevices", "bmp", "colorRampPalette", "dev.off", "jpeg",
                 "pdf", "png", "postscript")
      importFrom("graphics", "axis", "hist", "mtext", "par", "points",
                 "text")
      importFrom("stats", "ave", "dist", "hclust", "median", "quantile")
      importFrom("utils", "data", "read.table")
    to your NAMESPACE file.
    ```

