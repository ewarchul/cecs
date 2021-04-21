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
    stop("Missing argument: 'func_index' has to be provided !")
  }

  if (missing(x)) {
    stop("Missing argument: 'x' has to be provided !")
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
      problem = as.integer(func_index),
      input = as.double(x),
      row = as.integer(row),
      col = as.integer(col),
      output = double(row),
      PACKAGE = "cecs"
    )$output)
  } else {
    stop(
      stringr::str_glue(
        "Invalid argument: function index should be an integer between\\
        1 and 10!"
      )
    )
  }
}

##' CEC2019 interface
#'
#' @description
#' The R interface for CEC2019 100-Digit Challenge
#' Constrained Numerical Optimization benchmark.
#' Available dimensions are following: functions F1-F3 are available only for
#' (respective) dimensions 9, 16, and 18. Functions F4-F10 are available for
#' 10 dimensions.

#' @param func_index numeric index of optimisation problem from set set 1:10
#' @param x vector of numeric inputs for objective function
#' @return value of objective function for given input
#' @export
#' @useDynLib cecs

cec2019 <- function(func_index, x) {
  if (missing(func_index)) {
    stop("Missing argument: 'func_index' has to be provided !")
  }

  if (missing(x)) {
    stop("Missing argument: 'x' has to be provided !") 
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
    if (func_index == 1 && col != 9) {
      stop(
        stringr::str_glue(
          "Invalid argument: Function 1 is available only for 9 dimensions!"
        )
      )
    }
    else if (func_index == 2 && col != 16) {
      stop(
        stringr::str_glue(
          "Invalid argument: Function 2 is available only for 16 dimensions!"
        )
      )
    }
    else if (func_index == 3 && col != 18) {
      stop(
        stringr::str_glue(
          "Invalid argument: Function 3 is available only for 18 dimensions!"
        )
      )
    } 
    if ((func_index %in% 4:10) && col != 10) {
      stop(
        stringr::str_glue(
          "Invalid argument: Functions 4-10 are available only\\
          for 10 dimensions!"
        )
      )
    }
    extdatadir <- system.file("extdata/cec2019/", package = "cecs")
    if (extdatadir == "") {
      extdatadir <-
        unzip_data(download_data("cec2019"))
    }
    return(.C(
      "cecs",
      extdatadir = as.character(extdatadir),
      suite = as.character(""),
      cec = as.integer(19),
      problem = as.integer(func_index),
      input = as.double(x),
      row = as.integer(row),
      col = as.integer(col),
      output = double(row),
      PACKAGE = "cecs"
    )$output)
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
    stop("Missing argument: 'func_index' has to be provided !")
  }

  if (missing(x)) {
    stop("Missing argument: 'x' has to be provided !")
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
      problem = as.integer(func_index),
      input = as.double(x),
      row = as.integer(row),
      col = as.integer(col),
      output = double(row),
      PACKAGE = "cecs"
    )$output)
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
      problem = as.integer(func_index),
      input = as.double(x),
      row = as.integer(row),
      col = as.integer(col),
      output = double(row),
      PACKAGE = "cecs"
    )$output)
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
    stop("Missing argument: 'func_index' has to be provided !")
  }

  if (missing(x)) {
    stop("Missing argument: 'x' has to be provided !")
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
      problem = as.integer(func_index),
      input = as.double(x),
      row = as.integer(row),
      col = as.integer(col),
      output = double(row),
      PACKAGE = "cecs"
    )$output)
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
  cec2013::cec2013(func_index, x)
}
