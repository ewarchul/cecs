#' @export

cec2021 = function(i, x, suite) {
  if (base::missing(i)) base::stop("Missing argument; 'i' has to be provided !")

  if (base::missing(x)) base::stop("Missing argument; 'x' has to be provided !")

  if (base::is.numeric(i) && i >= 1 && i <= 30) {
    if (base::is.vector(x)) {
      row = 1
      col = base::length(x)
    } else if (base::is.matrix(x)) {
      row = base::nrow(x)
      col = base::ncol(x)
    } else {
      base::stop("x should be a vector or a matrix")
    } 

    if (!(col %in% c(2, 10, 20))) {
      base::stop("Invalid argument: Only 2, 10, and 20  dimensions/variables are allowed !")
    }
    extdatadir = base::system.file("extdata", package = "cec2021")
    f = base::.C("cec2021",
      extdatadir = as.character(extdatadir),
      i = base::as.integer(i), x = base::as.double(x), row = base::as.integer(row),
      col = base::as.integer(col), f = double(row), suite = base::as.character(suite),
      PACKAGE = "cec2021"
    )$f
  } else {
    base::stop("Invalid argument: 'i' should be an integer between 1 and 30 !")
  }
  return(f)
}