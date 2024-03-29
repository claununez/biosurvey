% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preselected_dist_mask.R
\name{preselected_dist_mask}
\alias{preselected_dist_mask}
\title{Helper to create objects to detect points to close to preselected sites}
\usage{
preselected_dist_mask(master, expected_points, space, variable_1 = NULL,
                      variable_2 = NULL, use_blocks = FALSE,
                      verbose = TRUE)
}
\arguments{
\item{master}{master_matrix object derived from function
\code{\link{prepare_master_matrix}} or master_selection object derived
from functions \code{\link{random_selection}},
\code{\link{uniformE_selection}}, or \code{\link{EG_selection}}.}

\item{expected_points}{(numeric) number of survey points (sites) to be
selected.}

\item{space}{(character) space in which the thinning will be performed. There
are two options available: "G", if it will be in geographic space, and
"E", if it will be in environmental space.}

\item{variable_1}{(character or numeric) name or position of the first
variable (x-axis). Default = NULL.}

\item{variable_2}{(character or numeric) name or position of the second
variable (y-axis). Default = NULL.}

\item{use_blocks}{(logical) whether or not to use block centroids instead of
all points when \code{space} = "E". Default = FALSE.}

\item{verbose}{(logical) whether or not to print messages about the process.
Default = TRUE.}
}
\value{
A list of two elements: the distance used to obtain \code{expected_points}
and a SpatVector object created from preselected_sites in
master.
}
\description{
Helper to create objects to detect points to close to preselected sites
}
\examples{
# Data
m_matrix_pre <- read_master(system.file("extdata/m_matrix_pre.rds",
                                        package = "biosurvey"))

# Running
pdm <- preselected_dist_mask(master = m_matrix_pre, expected_points = 20,
                             space = "G")
}
