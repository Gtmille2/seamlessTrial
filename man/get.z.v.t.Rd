% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.z.v.t.R
\name{get.z.v.t}
\alias{get.z.v.t}
\title{get.z.v.t function}
\usage{
get.z.v.t(data, selected, n1, N1, N)
}
\arguments{
\item{data}{data set with primary and secondary endpoint, and treatment}

\item{selected}{The selected treatment from the first look}

\item{n1}{Number of patients with primary endpoint available at first analysis}

\item{N1}{Number of patients with secondary endpoint available at first analysis}

\item{N}{The total number of patients in the trial}
}
\description{
gets z and v for all treatments at first look and for selected treatment at second look selection is prespecified or is of best treatment from first look if selected = 0
}
