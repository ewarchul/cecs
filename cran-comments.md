# Resubmission

This is a resubmission. In this verion I have:

* fixed memory leaks in original C code
* refactored original C code.

## Valgrind log

```
==646303== Memcheck, a memory error detector
==646303== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
==646303== Using Valgrind-3.15.0 and LibVEX; rerun with -h for copyright info
==646303== Command: /usr/lib/R/bin/exec/R -f testthat.R --restore --save --no-readline --vanilla
==646303== 

R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> library(testthat)
> library(cecs)
> 
> test_check("cecs")
trying URL 'http://home.elka.pw.edu.pl/~ewarchul/cec2014.zip'
Content type 'application/zip' length 2981864 bytes (2.8 MB)
==================================================
downloaded 2.8 MB

trying URL 'http://home.elka.pw.edu.pl/~ewarchul/cec2015.zip'
Content type 'application/zip' length 2034138 bytes (1.9 MB)
==================================================
downloaded 1.9 MB

trying URL 'http://home.elka.pw.edu.pl/~ewarchul/cec2017.zip'
Content type 'application/zip' length 3757812 bytes (3.6 MB)
==================================================
downloaded 3.6 MB

trying URL 'http://home.elka.pw.edu.pl/~ewarchul/cec2019.zip'
Content type 'application/zip' length 22326 bytes (21 KB)
==================================================
downloaded 21 KB

trying URL 'http://home.elka.pw.edu.pl/~ewarchul/cec2021.zip'
Content type 'application/zip' length 113654 bytes (110 KB)
==================================================
downloaded 110 KB

[ FAIL 0 | WARN 0 | SKIP 0 | PASS 836 ]
> 
> proc.time()
   user  system elapsed 
375.728   1.476 381.354 
==646303== 
==646303== HEAP SUMMARY:
==646303==     in use at exit: 89,629,364 bytes in 16,984 blocks
==646303==   total heap usage: 739,667 allocs, 722,683 frees, 913,745,917 bytes allocated
==646303== 
==646303== LEAK SUMMARY:
==646303==    definitely lost: 0 bytes in 0 blocks
==646303==    indirectly lost: 0 bytes in 0 blocks
==646303==      possibly lost: 0 bytes in 0 blocks
==646303==    still reachable: 89,629,364 bytes in 16,984 blocks
==646303==         suppressed: 0 bytes in 0 blocks
==646303== Reachable blocks (those to which a pointer was found) are not shown.
==646303== To see them, rerun with: --leak-check=full --show-leak-kinds=all
==646303== 
==646303== For lists of detected and suppressed errors, rerun with: -s
==646303== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

```

## Test environments
* local R installation, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 note


# Previous resubmission

This is a resubmission. In this verion I have:

* deleted github LICENSE file
* deleted "+ LICENSE" from DESCRIPTION 
* add trailing slash to URL.


## Test environments
* local R installation, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 3.6.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 note

# Previous resubmission

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
