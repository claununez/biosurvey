% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection_helpers.R
\name{point_thinning}
\alias{point_thinning}
\title{Helps in thinning points either in geographic or environmental space}
\usage{
point_thinning(data, x_column, y_column, thinning_distance, space,
               max_n_samples = 1, replicates = 10, set_seed = 1)
}
\arguments{
\item{data}{matrix or data.frame that contains at least two columns.}

\item{x_column}{(character) the name of the x-axis.}

\item{y_column}{(character) the name of the y-axis.}

\item{thinning_distance}{(numeric) distance for thinning. Units must be
selected according to the space, kilometers (km) for geographic and
Euclidean distances for environmental space.}

\item{space}{(character) space in which the thinning will be performed.
There are two options available: "G", if it will be in geographic space, and
"E", if it will be in environmental space.}

\item{max_n_samples}{(numeric) maximum number of samples to chose with most
points included. Default = 1.}

\item{replicates}{(numeric) number of thinning replicates. Default = 10.}

\item{set_seed}{(numeric) integer value to specify a initial seed.
Default = 1.}
}
\value{
A list with one or more elements, depending on \code{max_n_samples}. Each
element is a data.frame containing points retained after thinning. All
elements are different in at least one of the selected points.
}
\description{
Point thinning based on user-defined distances in geographic or
environmental space.
}
\examples{
# Data
m_matrix <- read_master(system.file("extdata/m_matrix.rds",
                                    package = "biosurvey"))
data1 <- m_matrix$data_matrix

# Thinning the points
thin <- point_thinning(data1, x_column = "Longitude", y_column = "Latitude",
                       thinning_distance = 200, space = "G",
                       max_n_samples = 1, replicates = 5, set_seed = 1)
}
