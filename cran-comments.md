# CRAN submission aroma.affymetrix 2.14.0
on 2015-10-24

Updates related to R / CRAN:

* Explicitly importing core R functions

Thanks in advance


## Notes not sent to CRAN
The package has been verified using `R CMD check --as-cran` on:

* Platform x86_64-pc-linux-gnu (64-bit):
  - R version 3.1.1 (2014-07-10)
  - R version 3.1.3 (2015-03-09)
  - R version 3.2.0 (2015-04-16)
  - R version 3.2.2 (2015-08-14)
  - R version 3.2.2 Patched (2015-10-19 r69550)
  - R Under development (unstable) (2015-10-23 r69563)

* Platform: x86_64-apple-darwin13.4.0 (64-bit):
  - R version 3.2.2 Patched (2015-10-22 r69556)

* Platform x86_64-w64-mingw32/x64 (64-bit):
  - R version 3.1.3 (2015-03-09)
  - R version 3.2.2 (2015-08-14)
  - R version 3.2.2 Patched (2015-10-19 r69550)

It has also verified using the <http://win-builder.r-project.org/> service.

Moreover, the updates cause no issues for any of the following
7 reverse dependency on CRAN and Bioconductor, which have been
tested with `R CMD check --as-cran`: ACNE 0.8.0, MPAgenomics 1.1.2,
NSA 0.0.32, PECA 1.6.0, Repitools 1.16.0 and TIN 1.2.0.
