% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/growthfd.bgs.R
\name{growthfd.bgs.apvs}
\alias{growthfd.bgs.apvs}
\title{Find apvs on growth curves}
\usage{
growthfd.bgs.apvs(age, velocity)
}
\arguments{
\item{age}{Vector of ages}

\item{velocity}{Matrix of velocity curves}

\item{ids}{Vector of individuals' ids}

\item{limits}{List of limits}
}
\value{
Vector of apv values
}
\description{
This function finds ages of maximal growth velocity on the velocity curves.
}
