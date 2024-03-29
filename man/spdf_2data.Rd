% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pam_helpers.R
\name{spdf_2data}
\alias{spdf_2data}
\title{Creates a data.frame of species' references from SpatVector}
\usage{
spdf_2data(spdf_object, spdf_grid, parallel = FALSE, n_cores = NULL)
}
\arguments{
\item{spdf_object}{SpatVector representing species' geographic
distributions. The data.frame associated with the object must contain a
column named "Species" to distinguish among features.}

\item{spdf_grid}{SpatVector of geographic grid for the region of interest
(output of function \code{\link{grid_from_region}}).}

\item{parallel}{(logical) whether to perform analyses in parallel.
Default = FALSE.}

\item{n_cores}{(numeric) number of cores to be used when \code{parallel} =
TRUE. The default, NULL, uses available cores - 1.}
}
\value{
A data.frame of species' found in distinct positions (defined with
identifiers); includes two columns: "ID" and "Species".
}
\description{
Creates a data.frame of species' references that contains
identifiers of position and species name, using a SpatVector representing
multiple species as input.
}
\examples{
# Data
species_data <- terra::vect(system.file("extdata/species_data.gpkg",
                                        package = "biosurvey"))
mx <- terra::vect(system.file("extdata/mx.gpkg", package = "biosurvey"))

# GRID
grid_reg <- grid_from_region(region = mx, cell_size = 100)

# Species data from polygons
sp_data <- spdf_2data(spdf_object = species_data, spdf_grid = grid_reg)
summary(sp_data)
}
