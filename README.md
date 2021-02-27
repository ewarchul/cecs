[![R-CMD-check](https://github.com/ewarchul/cecs/workflows/R-CMD-check/badge.svg)](https://github.com/ewarchul/cecs/actions)
[![](https://www.r-pkg.org/badges/version/cecs?color=green)](https://cran.r-project.org/package=cecs)
# cecs

Common R interface for CEC benchmarks.

<!--ts-->
   * [Descritpion](#description)
   * [Requirements](#requirements)
   * [Installation](#installation)
      * [Linux and Windows](#linux-and-windows)
      * [macOS](#macos)
   * [Content](#content)
     * [CEC 2021](#cec-2021)
     * [CEC 2017](#cec-2017)
     * [CEC 2014](#cec-2014)
     * [CEC 2013](#cec-2013)
<!--te-->

## Description

This repository contains **R** package with an interface for benchmark functions from Congress on Evolutionary Computations competitions. 
The implementation is based on the existing CRAN package, i.e `{cec2013}`, but is extended with CEC2014 and CEC2021. The interface for CEC2017 is taken from the non-CRAN package `{cec2017}`. More information about authorship is written below. 

## Requirements

* **R** in version 3.6.3

* `{stringr}` package in version 1.4.0

* `gcc` compiler in version 9.3.0

Package was written and tested on machine with Ubuntu 20.04.1 LTS focal OS and x86_64 CPU architecture.

## Installation

### Linux and Windows

#### CRAN (temporary archived)

Type bellow command in R interpreter:

```r
install.packages("cecs")
```

#### GitHub

Type bellow command in R interpreter:

```r
library(devtools)
devtools::install_github("ewarchul/cecs")
```

### macOS

#### CRAN (temporary archived)

```r
install.packages("cecs")
```

#### GitHub

Type bellow command in R interpreter:

```r
library(devtools)
devtools::install_github("ewarchul/cecs")
```

## Content

All benchmark functions were implemented in **C** by Jane Jing Liang [https://orcid.org/0000-0003-0811-0223](https://orcid.org/0000-0003-0811-0223).

I rearranged the codebase, i.e. separated interface from source, for better readability but any major changes were not done.

Benchmarks specifactions and necessary numeric data are available here: [http://home.elka.pw.edu.pl/~ewarchul/](http://home.elka.pw.edu.pl/~ewarchul/).


The package downloads numeric data from given website. For futher details, see source code documentation in `man/` directory.

### CEC 2021 

#### Authorship

> Authors: Eryk Warchulski

> License: GPL (>=3)

> Source: herein

### CEC 2017

#### Authorship 

> Authors: Dariusz Jagodzinski, 

> License: GPL (>=3)

> Source: [http://staff.elka.pw.edu.pl/~djagodzi/programy.html](http://staff.elka.pw.edu.pl/~djagodzi/programy.html)

### CEC 2014

#### Authorship

> Authors: Eryk Warchulski 

> License: GPL (>=3)

> Source: herein

### CEC 2013 

#### Authorship

> Authors: Mauricio Zambrano-Bigiarini, Yasser Gonzalez-Fernandez

> License: GPL (>=3)

> Source: [https://github.com/hzambran/cec2013](https://github.com/hzambran/cec2013)
