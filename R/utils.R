#' Extract compressed data
#'
#' @description
#' Function extracts CEC's data from ZIP archive.
#'
#' @param archivepath path to CEC's data

unzip_data <- function(archivepath) {
  exdir <- stringr::str_extract(archivepath, ".*extdata")
  cec <- stringr::str_extract(archivepath, "cec.\\d+")
  utils::unzip(
    archivepath,
    exdir = exdir
  )
  file.remove(archivepath)
  return(stringr::str_glue("{exdir}/{cec}"))
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

#' Download CEC data
#'
#' @description
#' Function downloads numeric data forspecified CEC
#' benchmark.
#' For further details, see \url{http://home.elka.pw.edu.pl/~ewarchul/}
#' @param cec name of benchmark

download_data <- function(cec) {
  url <- stringr::str_glue("http://home.elka.pw.edu.pl/~ewarchul/{cec}.zip")
  destfile <- stringr::str_glue(
    '{system.file("extdata/", package = "cecs")}{cec}.zip'
  )
  utils::download.file(url, destfile = destfile)
  return(destfile)
}
