% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconcile.R
\name{reconcile_train}
\alias{reconcile_train}
\title{function to train reconciliation model}
\usage{
reconcile_train(basef, real, lambda = 1, optimized = TRUE, step_wise = TRUE)
}
\arguments{
\item{basef}{list containing base probabilistic forecasts of each series}

\item{real}{dhts object containing corresponding observations of basef}

\item{lambda}{penality}

\item{optimized}{boolean indicating if optimization}

\item{step_wise}{boolean indicating if reconciliation step wisely.}
}
\value{
sw_res object or rec_mat object
}
\description{
function to train reconciliation model
}
