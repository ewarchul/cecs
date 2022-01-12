#' Get path to CEC data
#'
#' @description
#' If temporary directory with CEC's data exists then
#' function returns a filepath to data of specified
#' CEC version. Otherwise, function creates a temporary
#' directory, downloads an archive with data, unzip it,
#' and returns mentioned datapath.
#'
#' @param cec name of benchmark

get_cec_dirpath <- function(cec) {
  os <- get_os_name()
  tmp_dirpath <- cecs_tmp_dirpath(os)
  cec_datapath <- paste0(tmp_dirpath, cec)
  if (file.exists(cec_datapath)) {
    return(cec_datapath)
  } else {
    if (!dir.exists(tmp_dirpath)) {
      dir.create(tmp_dirpath)
    }
    download_data(cec, os)
    unzip_data(stringr::str_glue("{tmp_dirpath}{cec}.zip"), cec, os)
    return(cec_datapath)
  }
}

#' Remove CEC data
#'
#' @description
#' Function deletes directory with CEC's data.
#' Some text files are quite large (even 10MB)
#' and if one won't use CEC's data in the near
#' future this function allows to free the disk space.
#' @export

clean <- function() {
  os <- get_os_name()
  tmp_dirpath <- cecs_tmp_dirpath(os)
  if (base::dir.exists(tmp_dirpath)) {
    base::unlink(x = tmp_dirpath, recursive = TRUE)
  }
}

#' OS name
#'
#' @description
#' Helper function which returns OS name.

get_os_name <- function() {
  if (.Platform$OS.type == "windows") {
    "win"
  } else if (Sys.info()["sysname"] == "Darwin") {
    "osx"
  } else {
    "unix"
  }
}

#' Windows temporary directory name
#'
#' @description
#' Function returns a temporary directory name on Windows.

win_tmp_dirpath <- function() {
  return(stringr::str_extract(base::tempdir(), ".*Temp"))
}

#' cecs temporary directory name
#'
#' @description
#' Function returns a temporary directory name which
#' the package utilizes as a data storage on given
#' platform.
#' @param os_name OS name

cecs_tmp_dirpath <- function(os_name) {
  if (os_name == "win") {
    stringr::str_glue("{win_tmp_dirpath()}\\R-cecs-data\\")
  } else {
    stringr::str_glue("/tmp/R-cecs-data/")
  }
}

#' Download destination
#'
#' @description
#' Function returns filepath to the download
#' destination on given platform.
#' @param os_name OS name

destination_file <- function(os_name) {
  if (os_name == "win") {
    stringr::str_glue("{cecs_tmp_dirpath(os_name)}\\")
  } else {
    stringr::str_glue("{cecs_tmp_dirpath(os_name)}/")
  }
}

#' Download CEC data
#'
#' @description
#' Function downloads numeric data for specified CEC
#' benchmark.
#' For further details, see \url{https://github.com/ewarchul/cec/tree/main/data}
#' @param cec name of benchmark
#' @param os_name OS name

download_data <- function(cec, os_name) {
  url <- stringr::str_glue(
    "https://github.com/ewarchul/cec/raw/main/data/{cec}.zip"
  )
  dst <- destination_file(os_name)
  utils::download.file(
    url,
    destfile = stringr::str_glue("{dst}{cec}.zip")
  )
}

#' Extract compressed data
#'
#' @description
#' Function extracts CEC's data from ZIP archive.
#'
#' @param archivepath path to CEC's data
#' @param cec name of benchmark
#' @param os_name OS name

unzip_data <- function(archivepath, cec, os_name) {
  tmp_dirpath <- cecs_tmp_dirpath(os_name)
  utils::unzip(zipfile = archivepath, exdir = tmp_dirpath)
}
