# ACNE

Version: 0.8.1

## In both

*   checking CRAN incoming feasibility ... WARNING
    ```
    Maintainer: ‘Henrik Bengtsson <henrikb@braju.com>’
    
    Insufficient package version (submitted: 0.8.1, existing: 0.8.1)
    
    The Date field is over a month old.
    
    This build time stamp is over a month old.
    ```

# NSA

Version: 0.0.32

## In both

*   checking CRAN incoming feasibility ... WARNING
    ```
    Maintainer: ‘Maria Ortiz-Estevez <mortizest@gmail.com>’
    
    Insufficient package version (submitted: 0.0.32, existing: 0.0.32)
    
    The Title field should be in title case, current version then in title case:
    ‘Post-normalization of total copy numbers’
    ‘Post-Normalization of Total Copy Numbers’
    
    The Date field is over a month old.
    
    This build time stamp is over a month old.
    ```

*   checking package dependencies ... NOTE
    ```
    Depends: includes the non-default packages:
      ‘R.methodsS3’ ‘MASS’ ‘matrixStats’ ‘R.oo’ ‘R.utils’ ‘aroma.core’
      ‘aroma.affymetrix’ ‘DNAcopy’
    Adding so many packages to the search path is excessive and importing
    selectively is preferable.
    ```

*   checking top-level files ... NOTE
    ```
    Non-standard file/directory found at top level:
      ‘incl’
    ```

*   checking dependencies in R code ... NOTE
    ```
    Packages in Depends field not imported from:
      ‘DNAcopy’ ‘MASS’ ‘R.methodsS3’ ‘R.oo’ ‘aroma.affymetrix’ ‘aroma.core’
      ‘matrixStats’
      These packages need to be imported from (in the NAMESPACE file)
      for when this namespace is loaded but not attached.
    ```

*   checking S3 generic/method consistency ... NOTE
    ```
    Found the following apparent S3 methods exported but not registered:
      NSAByTotalAndFracB.matrix allocateOutputDataSets.NSANormalization
      allocateOutputDataSets.SNPsNormalization
      allocateOutputDataSets.SampleNormalization
      as.character.NSANormalization as.character.SNPsNormalization
      as.character.SampleNormalization findArraysTodo.NSANormalization
      findArraysTodo.SampleNormalization findUnitsTodo.SNPsNormalization
      fitNSA.matrix fitNSAcnPs.matrix getDataSets.NSANormalization
      getDataSets.SNPsNormalization getDataSets.SampleNormalization
      getName.NSANormalization getName.SNPsNormalization
      getName.SampleNormalization getOutputDataSets.NSANormalization
      getOutputDataSets.SNPsNormalization
      getOutputDataSets.SampleNormalization getPath.NSANormalization
      getPath.SNPsNormalization getPath.SampleNormalization
      getRootPath.NSANormalization getRootPath.SNPsNormalization
      getRootPath.SampleNormalization process.NSANormalization
      process.SNPsNormalization process.SampleNormalization
      sampleNByTotalAndFracB.numeric snpsNByTotalAndFracB.matrix
    See section ‘Registering S3 methods’ in the ‘Writing R Extensions’
    manual.
    ```

*   checking R code for possible problems ... NOTE
    ```
    ...
      ‘str’
      (/home/hb/repositories/aroma.affymetrix/revdep/checks/NSA/new/NSA.Rcheck/00_pkg_src/NSA/R/snpsNByTotalAndFracB.R:49)
    snpsNByTotalAndFracB.matrix: no visible global function definition for
      ‘rowAlls’
      (/home/hb/repositories/aroma.affymetrix/revdep/checks/NSA/new/NSA.Rcheck/00_pkg_src/NSA/R/snpsNByTotalAndFracB.R:54)
    snpsNByTotalAndFracB.matrix: no visible global function definition for
      ‘str’
      (/home/hb/repositories/aroma.affymetrix/revdep/checks/NSA/new/NSA.Rcheck/00_pkg_src/NSA/R/snpsNByTotalAndFracB.R:73)
    Undefined global functions or variables:
      AffymetrixCdfFile CNA Object approxfun aromaSettings byPath extend
      extractMatrix findUnitsTodo getAsteriskTags getChipType getFile
      getFullName getFullNames getGenomeInformation getName getNames
      getPath getPathname getPathnames getPositions getRam getRootPath
      getTags getUnitsOnChromosome hist median nbrOfFiles newInstance
      packageDescription rowAlls rowMedians segment setTags str throw trim
      verbose
    Consider adding
      importFrom("graphics", "hist")
      importFrom("stats", "approxfun", "median")
      importFrom("utils", "packageDescription", "str")
    to your NAMESPACE file.
    ```

# PECA

Version: 1.14.0

## In both

*   checking CRAN incoming feasibility ... NOTE
    ```
    Maintainer: ‘Tomi Suomi <tomi.suomi@utu.fi>’
    
    Package duplicated from https://bioconductor.org/packages/3.6/bioc
    
    The Title field should be in title case, current version then in title case:
    ‘Probe-level Expression Change Averaging’
    ‘Probe-Level Expression Change Averaging’
    
    This build time stamp is over a month old.
    ```

# REIDS

Version: 0.0.1

## In both

*   checking CRAN incoming feasibility ... WARNING
    ```
    Maintainer: ‘Marijke Van Moerbeke <marijke.vanmoerbeke@uhasselt.be>’
    
    Insufficient package version (submitted: 0.0.1, existing: 0.0.1)
    
    The Date field is over a month old.
    
    This build time stamp is over a month old.
    ```

# Repitools

Version: 1.24.0

## In both

*   checking CRAN incoming feasibility ... NOTE
    ```
    Maintainer: ‘Mark Robinson <mark.robinson@imls.uzh.ch>’
    
    Package duplicated from https://bioconductor.org/packages/3.6/bioc
    
    Found the following (possibly) invalid URLs:
      URL: http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/
        From: man/genQC.Rd
        Status: 403
        Message: Forbidden
    
    The Title field should be in title case, current version then in title case:
    ‘Epigenomic tools’
    ‘Epigenomic Tools’
    
    The Date field is not in ISO 8601 yyyy-mm-dd format.
    
    This build time stamp is over a month old.
    ```

# TIN

Version: 1.10.0

## In both

*   checking CRAN incoming feasibility ... NOTE
    ```
    Maintainer: ‘Bjarne Johannessen <bjajoh@rr-research.no>’
    
    Package duplicated from https://bioconductor.org/packages/3.6/bioc
    
    The Title field should be in title case, current version then in title case:
    ‘Transcriptome instability analysis’
    ‘Transcriptome Instability Analysis’
    
    The Date field is over a month old.
    
    This build time stamp is over a month old.
    ```

*   checking top-level files ... NOTE
    ```
    Non-standard file/directory found at top level:
      ‘doc’
    ```

