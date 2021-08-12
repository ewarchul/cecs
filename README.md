[![R-CMD-check](https://github.com/ewarchul/cecs/workflows/R-CMD-check/badge.svg)](https://github.com/ewarchul/cecs/actions)
[![](https://www.r-pkg.org/badges/version/cecs?color=green)](https://cran.r-project.org/package=cecs)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# cecs

Common R interface for CEC benchmarks.

<!--ts-->
   * [Description](#description)
   * [Requirements](#requirements)
   * [Installation](#installation)
      * [Linux and Windows](#linux-and-windows)
      * [macOS](#macos)
   * [Content](#content)
     * [CEC 2021](#cec-2021)
     * [CEC 2019](#cec-2019)
     * [CEC 2017](#cec-2017)
     * [CEC 2015](#cec-2015)
     * [CEC 2014](#cec-2014)
     * [CEC 2013](#cec-2013)
<!--te-->

## Description

This repository contains **R** package with an interface for benchmark functions from Congress on Evolutionary Computations competitions.
The implementation is based on the existing CRAN package, i.e `{cec2013}`, but is extended with CEC2014, CEC2015, CEC2019, and CEC2021.
The interface for CEC2017 is taken from the non-CRAN package `{cec2017}`.
More information about the authorship is written in `DESCRIPTION`. 

## Installation

### CRAN 

Type bellow command in R interpreter:

```r
install.packages("cecs")
```

### GitHub

Type bellow command in R interpreter:

```r
devtools::install_github("ewarchul/cecs")
```

## Content

All benchmark functions were implemented in **C** by Jane Jing Liang [https://orcid.org/0000-0003-0811-0223](https://orcid.org/0000-0003-0811-0223).

I rearranged the codebase, i.e. separated interface from source, changed global state management, etc.

Benchmarks specifications and necessary numeric data are available here: [http://home.elka.pw.edu.pl/~ewarchul/](http://home.elka.pw.edu.pl/~ewarchul/).

The package downloads numeric data from my website. For further details, see the source code documentation in `man/` directory.
