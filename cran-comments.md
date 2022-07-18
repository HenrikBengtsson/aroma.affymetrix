# CRAN submission aroma.affymetrix 3.2.1

on 2022-07-19

This submission addresses Rd issues detected by R-devel.

I've verified that this submission causes no issues for any of the 4 reverse package dependencies available on CRAN and Bioconductor.

Thanks in advance


## Notes not sent to CRAN

### R CMD check validation

The package has been verified using `R CMD check --as-cran` on:

| R version     | GitHub | R-hub    | mac/win-builder |
| ------------- | ------ | -------- | --------------- |
| 4.0.x         | L      |          |                 |
| 4.1.x         | L M W  |          |                 |
| 4.2.x         | L M W  | . M M1 . |  . W            |
| devel         | L M W  | .        |    W            |

*Legend: OS: L = Linux, M = macOS, M1 = macOS M1, W = Windows*


R-hub checks:

```r
res <- rhub::check(platform = c(
#  "debian-clang-devel", "debian-gcc-patched", "linux-x86_64-centos-epel",
  "macos-highsierra-release-cran", "macos-m1-bigsur-release"
#  "windows-x86_64-release"
))
print(res)
```

gives

```
── aroma.affymetrix 3.2.1: WARNING

  Build ID:   aroma.affymetrix_3.2.1.tar.gz-f5c48e38f1a04b5a87a304f4030d8d1e
  Platform:   macOS 10.13.6 High Sierra, R-release, CRAN's setup
  Submitted:  10m 24.5s ago
  Build time: 10m 22.3s

❯ checking whether package ‘aroma.affymetrix’ can be installed ... WARNING
  See below...

❯ checking installed package size ... NOTE
    installed size is  6.1Mb
    sub-directories of 1Mb or more:
      help          1.1Mb
      R             3.1Mb
      testScripts   1.1Mb

0 errors ✔ | 1 warning ✖ | 1 note ✖

── aroma.affymetrix 3.2.1: NOTE

  Build ID:   aroma.affymetrix_3.2.1.tar.gz-9da1b1f8dadf47d3ab61fe6f3de1099a
  Platform:   Apple Silicon (M1), macOS 11.6 Big Sur, R-release
  Submitted:  10m 24.5s ago
  Build time: 8m 57s

❯ checking installed package size ... NOTE
    installed size is  6.1Mb
    sub-directories of 1Mb or more:
      help          1.1Mb
      R             3.0Mb
      testScripts   1.1Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```
