% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cec.R
\name{cec2022}
\alias{cec2022}
\title{CEC2022 interface}
\usage{
cec2022(func_index, x)
}
\arguments{
\item{func_index}{numeric index of optimisation problem from set 1:12}

\item{x}{vector of numeric inputs for objective function}
}
\value{
value of objective function for given input
}
\description{
The R interface for CEC2022 Single Objective Bound
Constrained Numerical Optimization benchmark.
Available dimensions are following: (10, 20).
}
