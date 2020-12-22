# File cec2017.R
# Copyright 2019 Dariusz Jagodziński <d.jagodzinski@elka.pw.edu.pl>
# Based on: https://CRAN.R-project.org/package=cec2013
# Distributed under GPL 3 or later

################################################################################
#                                 CEC 2017                                     #
################################################################################
# Purpose   : Evaluate a CEC-2017 benchamark function on a user-defined para-  #
#             meter set                                                        #
################################################################################
# i: integer in [1, 30], with the number of the CEC2017 benchmark function to  #
#    be evalauated                                                             #
# x: numeric, with the parameter set to be evaluated in the benchmark function #
#    Its length MUST be in [2, 10, 20, 30, 50, 100]                            #
################################################################################

cec2017 <- function (i, x) {

  if (missing(i)) stop("Missing argument; 'i' has to be provided !")

  if (missing(x)) stop("Missing argument; 'x' has to be provided !")

  if (is.numeric(i) && i >= 1 && i <= 30) {
    if (is.vector(x)) {
      row <- 1; col <- length(x)
    } else if (is.matrix(x)) {
      row <- nrow(x); col <- ncol(x)
    } else {
      stop("x should be a vector or a matrix")
    } # ELSE end

    if (!(col %in% c(10, 20, 30, 50, 100))) {
      stop("Invalid argument: Only 10, 20, 30, 50 and 100 dimensions/variables are allowed !")
    }
    extdatadir <- system.file("extdata", package = "cec2017")
    f <- .C("cec2017", extdatadir = as.character(extdatadir),
            i = as.integer(i), x = as.double(x), row = as.integer(row),
            col = as.integer(col), f = double(row),
            PACKAGE = "cec2017")$f
  } else stop("Invalid argument: 'i' should be an integer between 1 and 30 !")

  return(f)
}
