% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/growthModelFunctions.R
\name{crossValidate}
\alias{crossValidate}
\title{Turn observedLength values to NA for a percentage of the observations}
\usage{
crossValidate(d, runCrossValidationTF)
}
\arguments{
\item{d}{a dataframe}

\item{runCrossValidation}{boolean for running cross validation or not}
}
\value{
a data frame with observedLength set to NA for percentLeftOut observations
}
\description{
Turn observedLength values to NA for a percentage of the observations
}
