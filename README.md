# cecs

Common R interface for CEC benchmarks.

## Description

This repository contains `R` package with interface for benchmark functions from Congress on Evolutionary Computations competitions. 
The implementation is based on existing CRAN package, i.e {cec2013}, but extends it with CEC2014 and CEC2021. The interface for CEC2017 is taken from non-CRAN package {cec2017}. More information about authorship is written below. 

### Requirements

* `R` in version 3.6.3

* `{stringr}` package in version 1.4.0

* `gcc` compiler in version 9.3.0

Package was written and tested on machine with Ubuntu 20.04.1 LTS focal OS and x86_64 CPU architecture.

## Installation

The easiest way to install {cecs} package is to use {devtools} package and run following command from `R` CLI:

```R
devtools::install_github("ewarchul/cecs")
```

After the installation run unit tests to check if the package works correctly:

```R
devtools::test()
```

## Content

All benchmark functions were implemented in `C` by Jane Jing Liang [https://orcid.org/0000-0003-0811-0223](https://orcid.org/0000-0003-0811-0223).

I rearranged the codebase, i.e. separated interface from source, for better readability but any major changes were not done.

Benchmarks specifactions are included in `inst/` directory as a PDF files.

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

## Todos

Due to low quality of implementation, CEC2005 is not included in the package.

It has to be rewritten and if you have enough time to do it -- feel free to contribute.





