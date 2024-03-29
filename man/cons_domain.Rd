% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dhts.R
\name{cons_domain}
\alias{cons_domain}
\title{construct domain}
\usage{
cons_domain(domain_bts, s_mat, coherent = TRUE, node_names = NULL)
}
\arguments{
\item{domain_bts}{domain matrix used in creating dhts object.}

\item{s_mat}{summing matrix}

\item{coherent}{logical indicating if produced domain is coherent or not.}
}
\value{
domain object
}
\description{
Method for constructing coherent domain vector given domain matrix of all 
bottom-level series and the summing matrix.
}
