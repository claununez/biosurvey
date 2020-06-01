% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_documentation.R
\docType{data}
\name{m_matrix}
\alias{m_matrix}
\title{Example of object obtained from using the master_matrix function}
\format{A list of 3 elements:
\describe{
\item{m_matrix}{data frame wiht 6276 rows and 10 columns}
\item{polygon}{object of class SpatialPolygonsDataFrame}
\item{PCA_results}{list of length 5}
}}
\usage{
m_matrix
}
\description{
A list of 3 elements (m_matrix, polygon, PCA_results). See function
\code{\link{master_matrix}}.
}
\examples{
data("m_matrix", package = "biosurvey")

summary(m_matrix)
}
\keyword{datasets}