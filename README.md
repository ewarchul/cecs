[![R-CMD-check](https://github.com/ewarchul/cecs/workflows/R-CMD-check/badge.svg)](https://github.com/ewarchul/cecs/actions)
[![](https://www.r-pkg.org/badges/version/cecs?color)](https://cran.r-project.org/package=cecs)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# cecs

Common R interface for CEC benchmarks.

<!--ts-->
   * [Description](#description)
   * [Requirements](#requirements)
   * [Installation](#installation)
   * [Content](#content)
<!--te-->

## Description

This repository contains **R** package with an interface for benchmark functions from Congress on Evolutionary Computations competitions.
The implementation is based on the existing CRAN package, i.e `{cec2013}`, but is extended with CEC2014, CEC2015, CEC2019, and CEC2021.
The interface for CEC2017 is taken from the non-CRAN package `{cec2017}`.
More information about the authorship is written in `DESCRIPTION`. 

## Installation

### CRAN 

Type bellow command in the **R** interpreter:

```r
install.packages("cecs")
```

### GitHub

Type bellow command in the **R** interpreter:

```r
require(devtools)
devtools::install_github("ewarchul/cecs")
```

## Content

All benchmark functions were implemented in **C** by Jane Jing Liang [https://orcid.org/0000-0003-0811-0223](https://orcid.org/0000-0003-0811-0223).

I rearranged the codebase, i.e. separated interface from source, changed global state management, etc.

Benchmarks specifications and necessary numeric data are available here: [http://home.elka.pw.edu.pl/~ewarchul/](http://home.elka.pw.edu.pl/~ewarchul/).

The package downloads numeric data from my website. For further details, see the source code documentation in `man/` directory.
