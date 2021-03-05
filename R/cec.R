#' CEC2021 interface
#'
#' @description
#' The R interface for CEC2021 Single Objective Bound
#' Constrained Numerical Optimization benchmark.
#' Available dimensions are following: (10, 20).
#'
#' @param func_index numeric index of optimisation problem from set 1:10
#' @param x vector of numeric inputs for objective function
#' @param suite one of the suite in CEC2021 benchmark
#' (basic, shift, rot, bias, shift_rot, bias_rot, bias_shift, bias_shift_rot)
#' @return value of objective function for given input
#' @export
#' @useDynLib cecs

cec2021 <- function(func_index, x, suite) {
  if (missing(func_index)) {
    stop("Missing argument; 'func_index' has to be provided !")
  }

  if (missing(x)) {
    stop("Missing argument; 'x' has to be provided !")
  }
  if (is.numeric(func_index) && func_index >= 1 && func_index <= 10) {
    if (is.vector(x)) {
      row <- 1
      col <- length(x)
    } else if (is.matrix(x)) {
      row <- nrow(x)
      col <- ncol(x)
    } else {
      stop("x should be a vector or a matrix")
    }
    if (!(col %in% c(10, 20))) {
      stop(
        stringr::str_glue(
          "Invalid argument: Only 10, 20 dimensions/variables are allowed!"
        )
      )
    }
    if (!(suite %in% c(
      "basic",
      "shift",
      "rot",
      "bias",
      "shift_rot",
      "bias_rot",
      "bias_shift",
      "bias_shift_rot"
    ))) {
      stop(
        stringr::str_glue(
          "Invalid argument: Only 10, 20 dimensions/variables are allowed!"
        )
      )
    }
    extdatadir <- system.file("extdata/cec2021/", package = "cecs")
    if (extdatadir == "") {
      extdatadir <-
        unzip_data(download_data("cec2021"))
    }
    return(.C(
      "cecs",
      extdatadir = as.character(extdatadir),
      suite = as.character(suite),
      cec = as.integer(21),
      i = as.integer(func_index),
      x = as.double(x),
      row = as.integer(row),
      col = as.integer(col),
      f = double(row),
      PACKAGE = "cecs"
    )$f)
  } else {
    stop(
      stringr::str_glue(
        "Invalid argument: function index should be an integer between\\
        1 and 10!"
      )
    )
  }
}

##' CEC2017 interface
#'
#' @description
#' The R interface for CEC2017 Single Objective Bound
#' Constrained Numerical Optimization benchmark.
#' Available dimensions are following: (10, 30, 50, 100).
#'
#' @param func_index numeric index of optimisation problem from set 1:30
#' @param x vector of numeric inputs for objective function
#' @return value of objective function for given input
#' @source http://staff.elka.pw.edu.pl/~djagodzi/programy.html
#' @export
#' @useDynLib cecs

cec2017 <- function(func_index, x) {
  if (missing(func_index)) {
    stop("Missing argument; 'func_index' has to be provided !")
  }

  if (missing(x)) {
    stop("Missing argument; 'x' has to be provided !")
  }
  if (is.numeric(func_index) && func_index >= 1 && func_index <= 30) {
    if (is.vector(x)) {
      row <- 1
      col <- length(x)
    } else if (is.matrix(x)) {
      row <- nrow(x)
      col <- ncol(x)
    } else {
      stop("x should be a vector or a matrix")
    }
    if (!(col %in% c(10, 30, 50, 100))) {
      stop(
        stringr::str_glue(
          "Invalid argument: Only 10, 30, 50, 100\\
          dimensions/variables are allowed !"
        )
      )
    }
    extdatadir <- system.file("extdata/cec2017/", package = "cecs")
    if (extdatadir == "") {
      extdatadir <-
        unzip_data(download_data("cec2017"))
    }
    return(.C(
      "cecs",
      extdatadir = as.character(extdatadir),
      suite = as.character(""),
      cec = as.integer(17),
      i = as.integer(func_index),
      x = as.double(x),
      row = as.integer(row),
      col = as.integer(col),
      f = double(row),
      PACKAGE = "cecs"
    )$f)
  } else {
    stop(
      stringr::str_glue(
        "Invalid argument: function index should be an integer between\\
        1 and 30!"
      )
    )
  }
}

##' CEC2015 interface
#'
#' @description
#' The R interface for CEC2015 Single Objective Bound
#' Constrained Numerical Optimization benchmark.
#' Available dimensions are following: (10, 30, 50, 100)
#'
#' @param func_index numeric index of optimisation problem from set set 1:15
#' @param x vector of numeric inputs for objective function
#' @return value of objective function for given input
#' @export
#' @useDynLib cecs

cec2015 <- function(func_index, x) {
  if (missing(func_index)) {
    stop("Missing argument; 'func_index' has to be provided!")
  }

  if (missing(x)) {
    stop("Missing argument; 'x' has to be provided!")
  }
  if (is.numeric(func_index) && func_index >= 1 && func_index <= 15) {
    if (is.vector(x)) {
      row <- 1
      col <- length(x)
    } else if (is.matrix(x)) {
      row <- nrow(x)
      col <- ncol(x)
    } else {
      stop("x should be a vector or a matrix")
    }
    if (!(col %in% c(10, 30, 50, 100))) {
      stop(
        stringr::str_glue(
          "Invalid argument: Only 10, 30, 50, 100\\
          dimensions/variables are allowed!"
        )
      )
    }
    extdatadir <- system.file("extdata/cec2015/", package = "cecs")
    if (extdatadir == "") {
      extdatadir <-
        unzip_data(download_data("cec2015"))
    }
    return(.C(
      "cecs",
      extdatadir = as.character(extdatadir),
      suite = as.character(""),
      cec = as.integer(15),
      i = as.integer(func_index),
      x = as.double(x),
      row = as.integer(row),
      col = as.integer(col),
      f = double(row),
      PACKAGE = "cecs"
    )$f)
  } else {
    stop(
      stringr::str_glue(
        "Invalid argument: function index should be an integer between\\
        1 and 15!"
      )
    )
  }
}


##' CEC2014 interface
#'
#' @description
#' The R interface for CEC2014 Single Objective Bound
#' Constrained Numerical Optimization benchmark.
#' Available dimensions are following: (10, 20, 30, 50, 100).
#'
#' @param func_index numeric index of optimisation problem from set set 1:30
#' @param x vector of numeric inputs for objective function
#' @return value of objective function for given input
#' @export
#' @useDynLib cecs

cec2014 <- function(func_index, x) {
  if (missing(func_index)) {
    stop("Missing argument; 'func_index' has to be provided !")
  }

  if (missing(x)) {
    stop("Missing argument; 'x' has to be provided !")
  }
  if (is.numeric(func_index) && func_index >= 1 && func_index <= 30) {
    if (is.vector(x)) {
      row <- 1
      col <- length(x)
    } else if (is.matrix(x)) {
      row <- nrow(x)
      col <- ncol(x)
    } else {
      stop("x should be a vector or a matrix")
    }
    if (!(col %in% c(10, 20, 30, 50, 100))) {
      stop(
        stringr::str_glue(
          "Invalid argument: Only 10, 20, 30, 50, 100\\
          dimensions/variables are allowed !"
        )
      )
    }
    extdatadir <- system.file("extdata/cec2014/", package = "cecs")
    if (extdatadir == "") {
      extdatadir <-
        unzip_data(download_data("cec2014"))
    }
    return(.C(
      "cecs",
      extdatadir = as.character(extdatadir),
      suite = as.character(""),
      cec = as.integer(14),
      i = as.integer(func_index),
      x = as.double(x),
      row = as.integer(row),
      col = as.integer(col),
      f = double(row),
      PACKAGE = "cecs"
    )$f)
  } else {
    stop(
      stringr::str_glue(
        "Invalid argument: function index should be an integer between\\
        1 and 30!"
      )
    )
  }
}

#' CEC2013 interface
#'
#' @description
#' The R interface for CEC2013 Single Objective Bound
#' Constrained Numerical Optimization benchmark.
#' Available dimensions are following:
#' (2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100).
#'
#' @param func_index numeric index of optimisation problem from set 1:28
#' @param x vector of numeric inputs for objective function
#' @return value of objective function for given input
#' @source https://github.com/hzambran/cec2013
#' @export
#' @useDynLib cecs

cec2013 <- function(func_index, x) {
  if (missing(func_index)) {
    stop("Missing argument; 'func_index' has to be provided !")
  }

  if (missing(x)) {
    stop("Missing argument; 'x' has to be provided !")
  }
  if (is.numeric(func_index) && func_index >= 1 && func_index <= 28) {
    if (is.vector(x)) {
      row <- 1
      col <- length(x)
    } else if (is.matrix(x)) {
      row <- nrow(x)
      col <- ncol(x)
    } else {
      stop("x should be a vector or a matrix")
    }
    if (!(col %in% c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))) {
      stop(
        stringr::str_glue(
          "Invalid argument: Only 2, 5, seq(10, 100, 10)\\
          dimensions/variables are allowed !"
        )
      )
    }
    extdatadir <- system.file("extdata/cec2013/", package = "cecs")
    if (extdatadir == "") {
      extdatadir <-
        unzip_data(download_data("cec2013"))
    }
    return(.C(
      "cecs",
      extdatadir = as.character(extdatadir),
      suite = as.character(""),
      cec = as.integer(13),
      i = as.integer(func_index),
      x = as.double(x),
      row = as.integer(row),
      col = as.integer(col),
      f = double(row),
      PACKAGE = "cecs"
    )$f)
  } else {
    stop(
      stringr::str_glue(
        "Invalid argument: function index should be an integer between\\
        1 and 28!"
      )
    )
  }
}
