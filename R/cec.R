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

cec2021 <- function(func_index, x, suite) {
  cec(
    func_index,
    x,
    cec = "cec2021",
    max_func_index = 10,
    dims = c(10, 20),
    suite = suite
  )
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

cec2017 <- function(func_index, x) {
  cec(
    func_index,
    x,
    cec = "cec2017",
    max_func_index = 30,
    dims = c(10, 30, 50, 100),
    suite = NULL
  )
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

cec2014 <- function(func_index, x) {
  cec(
    func_index,
    x,
    cec = "cec2014",
    max_func_index = 30,
    dims = c(10, 20, 30, 50, 100),
    suite = NULL
  )
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

cec2013 <- function(func_index, x) {
  cec(
    func_index,
    x,
    cec = "cec2013",
    max_func_index = 28,
    dims = c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
    suite = NULL
  )
}

#' CEC interface
#'
#' @description
#' The common interface for all available benchmark from CEC.
#'
#' @param func_index numeric index of optimisation problem in given benchmark
#' @param x vector of numeric inputs for objective function
#' @param cec name of benchmark
#' @param max_func_index biggest index of optimisation
#' problem in given benchmark
#' @param dims vector of available dimensionalities
#' @param suite one of the suite in CEC2021
#' benchmark (basic, shift, rot, bias,
#' shift_rot, bias_rot, bias_shift, bias_shift_rot)
#' @return value of objective function for given input
#' @useDynLib cecs

cec <- function(func_index, x, cec, max_func_index, dims, suite = NULL) {
  if (base::missing(func_index))
    base::stop("Missing argument; 'func_index' has to be provided !")

  if (base::missing(x))
    base::stop("Missing argument; 'x' has to be provided !")

  if (
    base::is.numeric(func_index) &&
      func_index >= 1 &&
      func_index <= max_func_index
  ) {
    if (base::is.vector(x)) {
      row <- 1
      col <- base::length(x)
    } else if (base::is.matrix(x)) {
      row <- base::nrow(x)
      col <- base::ncol(x)
    } else {
      base::stop("x should be a vector or a matrix")
    }
    if (!(col %in% dims)) {
      base::stop(
        stringr::str_interp(
          "Invalid argument: only ${dims} dimensions/variables are allowed !"
        )
      )
    }
    datapath <- base::paste0("extdata/", cec)
    extarchive <- system.file(
      base::paste0(datapath, ".zip"),
      package = "cecs"
    )
    unzip_data(extarchive, cec)
    extdatadir <- system.file(datapath, package = "cecs")
    f <- base::.C(
      "cecs",
      extdatadir = as.character(extdatadir),
      cec = base::as.character(cec),
      i = base::as.integer(func_index),
      x = base::as.double(x),
      row = base::as.integer(row),
      col = base::as.integer(col),
      f = double(row),
      suite = base::as.character(suite),
      PACKAGE = "cecs"
    )$f
  } else {
    base::stop(
      stringr::str_interp(
        "Invalid argument: 
        'func_index' should be an integer between 1 and ${max_func_index} !"
      )
    )
  }
  return(f)
}

#' Extract compressed data
#'
#' @description
#' Function extracts CEC's data from ZIP archive.
#' If data is already extracted function does nothing.
#'
#' @param datapath path to CEC's data
#' @param cec name of benchmark

unzip_data <- function(datapath, cec) {
  if (!base::dir.exists(stringr::str_sub(datapath, end = -5))) {
    utils::unzip(
      datapath,
      exdir = system.file("extdata/", package = "cecs")
    )
  }
}

#' Remove TXT data
#'
#' @description
#' Function deletes decompressed directory with
#' CEC's data. Some text files are quite large (even 10MB)
#' and if one won't use specific CEC version in
#' the near future this function allows to free the disk space.
#' @export

clean <- function() {
  cecs <- c(
    "cec2013",
    "cec2014",
    "cec2017",
    "cec2021"
  )
  purrr::walk(cecs, function(cec) {
    datadir <- system.file(
      paste0("extdata/", cec),
      package = "cecs"
    )
    if (dir.exists(datadir)) {
      base::unlink(datadir, recursive = TRUE)
    }
  })
}
