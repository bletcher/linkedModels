% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detectionModelFunctions.R
\name{getDensities}
\alias{getDensities}
\title{Get (estimate) densities by section, river, year based on num sammples and p from the density model}
\usage{
getDensities(dddd, dd, meanOrIter = "mean", sampleToUse = sampleToUse)
}
\arguments{
\item{dddd}{the input data frame for the detection model run}

\item{dd}{the output data frame for the detection model run}

\item{meanOrIter}{Whether model run output is the mean ('mean') of all iterations or a single iteration ('iter'). Default is 'mean'.}

\item{sampleToUse}{if meanOrIter == 'iter', which iteration to use. Assumes use of chain 1.}

\item{chainToUse}{if meanOrIter == 'iter', which chain to use. Default is chain 1.}
}
\value{
a data frame
}
\description{
Get (estimate) densities by section, river, year based on num sammples and p from the density model
}
