% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/short_helpers.R
\name{match_rformat}
\alias{match_rformat}
\title{Helper function to find raster extension}
\usage{
match_rformat(format)
}
\arguments{
\item{format}{(character) any of the format names allowed for raster objects.
Options are: "GTiff", "ascii", "EHdr", "SAGA", "IDRISI", "CDF", "ENVI",
and "HFA".}
}
\value{
Raster extension according to format type.
}
\description{
Helper function to find raster extension
}
\examples{
match_rformat("GTiff")
}
