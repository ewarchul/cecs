# Resubmission

This is a resubmission of previously archived package. The package was archived
due to the presence of the obsolete header <malloc.h> in C source code.

In this version I have:

* removed <malloc.h> header
* extened package with CEC 2015 and CEC 2019 
* improved performance.

## Test environments
* local R installation, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 note

# Previous resubmission

This is a resubmission. In this verion I have:

* moved data files from inst/extdata directory to the external source (http://home.elka.pw.edu.pl/~ewarchul/)
* provided functions to download and extract data 
* provided function to delete extracted archive data in inst/extdata directory 
* added reference to specifications of benchmark in the DESCRIPTION.

## Test environments
* local R installation, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

# Previous resubmission

This is a resubmission. In this verion I have:

* fixed the title field: reduce the length of the title, set title in title case.
* extended Autohrs@R field by copyright holders of source codes of cec2013 and cec2017


## Test environments
* local R installation, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
