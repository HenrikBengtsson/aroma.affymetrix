# CRAN submission aroma.affymetrix 3.1.1

on 2018-04-16

This submission addresses an R CMD check WARNING on R-devel due to
a 'break' statement in a non-loop construct.

Note: This submission does _not_ address the ERROR on R-devel related
to propotypical updates in S3-lookup methods.
This was approved by CRAN by private email.

I've verified that this submission causes no issues for any of the
6 reverse (non-recursive) package dependencies available on CRAN
and Bioconductor.

Thanks in advance.


## Notes not sent to CRAN

The package has been verified using `R CMD check --as-cran` on:

* Platform x86_64-apple-darwin13.4.0 (64-bit) [Travis CI]:
  - R version 3.3.3 (2017-01-27)
  - R version 3.4.4 (2018-03-15)
  
* Platform x86_64-unknown-linux-gnu (64-bit) [Travis CI]:
  - R version 3.3.3 (2017-01-27)
  - R version 3.4.4 (2017-01-27)
  - R Under development (unstable) (2018-04-16 r74608)

* Platform x86_64-pc-linux-gnu (64-bit):
  - R version 3.2.0 (2015-04-16)
  - R version 3.4.4 (2018-03-15)

* Platform x86_64-pc-linux-gnu (64-bit) [r-hub]:
  - R version 3.4.4 (2018-03-15)
  - R Under development (unstable) (2018-04-14 r74601)

* Platform i386-w64-mingw32 (32-bit) [Appveyor CI]:
  - R Under development (unstable) (2018-04-15 r74605)

* Platform x86_64-w64-mingw32/x64 (64-bit) [Appveyor CI]:
  - R version 3.4.4 (2018-03-15)
  - R Under development (unstable) (2018-04-15 r74605)

* Platform x86_64-w64-mingw32/x64 (64-bit) [win-builder]:
  - R version 3.4.4 (2018-03-15)
  - R version 3.5.0 RC (2018-04-15 r74605) [known, accepted ERROR]

