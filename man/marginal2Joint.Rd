% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{marginal2Joint}
\alias{marginal2Joint}
\title{convert marginal distributions into joint distribution assuming independence}
\usage{
marginal2Joint(x, hier, method = c("ind", "bu"))
}
\arguments{
\item{x}{list of distributions of all series}

\item{hier}{dhts object}
}
\value{
joint distribution matrix
}
\description{
convert marginal distributions into joint distribution assuming independence
}
