% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Classes.R
\name{PAM_CS}
\alias{PAM_CS}
\alias{new_PAM_CS}
\title{Constructor of S3 objects of class PAM_CS}
\usage{
new_PAM_CS(Species = NA, Sites_cells = NA, Beta_W = NA, Spearman_cor = NA,
           Theoretical_boundaries = list(x = NA, y = NA),
           Richness_normalized = NA, Dispersion_field_normalized = NA,
           S_significance_id = NA, Randomized_DF = matrix())
}
\arguments{
\item{Species}{(numeric) species name. Default = NA.}

\item{Sites_cells}{(numeric) number of sites or cells. Default = NA.}

\item{Beta_W}{(numeric) value of Whittaker's Beta. Default = NA.}

\item{Spearman_cor}{(numeric) value of Spearman's correlation. Default = NA.}

\item{Theoretical_boundaries}{list of theoretical boundaries for the values.
Default = NA.}

\item{Richness_normalized}{(numeric) values of normalized richness.
Default = NA.}

\item{Dispersion_field_normalized}{(numeric) values of normalized dispersion
field. Default = NA.}

\item{S_significance_id}{(numeric) values indicating statistical significance
of the normalized dispersion field. Default = NA.}

\item{Randomized_DF}{matrix of values resulted from randomizing matrices.
Default = NA.}
}
\value{
An object of class \code{PAM_CS}.
}
\description{
Constructor of S3 objects of class PAM_CS
}
