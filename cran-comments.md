# CRAN submission aroma.affymetrix 3.2.0

on 2019-06-22

This submission fixes long-standing R CMD check ERRORs that appear on the CRAN MS Windows servers.

I've verified that this submission causes no issues for any of the 5 reverse package dependencies available on CRAN and Bioconductor.

Thanks in advance


## Notes not sent to CRAN

### R CMD check --as-cran validation

The package has been verified using `R CMD check --as-cran` on:

* Platform x86_64-apple-darwin15.6.0 (64-bit) [Travis CI]:
  - R version 3.5.3 (2019-03-11)
  - R version 3.6.0 (2019-04-26)

* Platform x86_64-unknown-linux-gnu (64-bit) [Travis CI]:
  - R version 3.5.3 (2017-01-27) [sic!]
  - R version 3.6.0 (2017-01-27) [sic!]
  - R Under development (unstable) (2019-06-17 r76707)

* Platform x86_64-pc-linux-gnu (64-bit):
  - R version 3.2.0 (2015-04-16)
  - R version 3.3.0 (2016-05-03)
  - R version 3.4.0 (2017-04-21)
  - R version 3.5.0 (2018-04-23)
  - R version 3.6.0 (2019-04-26)
  - 3.6.0 Patched (2019-05-31 r76629)
  
* Platform x86_64-pc-linux-gnu (64-bit) [r-hub]:
  - R version 3.6.0 (2019-04-26)
  - R Under development (unstable) (2019-06-14 r76701)

* Platform i686-pc-linux-gnu (32-bit):
  - R version 3.4.4 (2018-03-15)

* Platform i386-pc-solaris2.10 (32-bit) [r-hub]:
  - R version 3.5.0 Patched (2018-04-30 r74674)

* Platform x86_64-w64-mingw32 (64-bit) [r-hub]:
  - R Under development (unstable) (2019-06-11 r76694)

* Platform i386-w64-mingw32 (32-bit) [Appveyor CI]:
  - R Under development (unstable) (2019-06-10 r76691)
  
* Platform x86_64-w64-mingw32/x64 (64-bit) [Appveyor CI]:
  - R Under development (unstable) (2019-06-10 r76691)

* Platform x86_64-w64-mingw32/x64 (64-bit) [win-builder]:
  - R version 3.6.0 (2019-04-26)
  - R Under development (unstable) (2019-06-20 r76728)

Not available:

* Platform x86_64-w64-mingw32/x64 (64-bit) [Appveyor CI]:
  - R version 3.6.0 (2019-04-26)
  REASON: Odd AppVeyor error on "Error: Bioconductor version not changed"
