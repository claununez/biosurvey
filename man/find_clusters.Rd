% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectionEG_helpers.R
\name{find_clusters}
\alias{find_clusters}
\title{Detection of clusters in 2D spaces}
\usage{
find_clusters(data, x_column, y_column, space,
              cluster_method = "hierarchical", n_k_means = NULL,
              split_distance = NULL)
}
\arguments{
\item{data}{matrix or data.frame that contains at least two columns.}

\item{x_column}{(character) the name of the x-axis.}

\item{y_column}{(character) the name of the y-axis.}

\item{space}{(character) space in which the thinning will be performed.
There are two options available: "G", if it will be in the geographic space,
and "E", if it will be in the environmental space.}

\item{cluster_method}{(character) name of the method to be used for detecting
clusters. Options are "hierarchical" and "k-means"; default = "hierarchical".}

\item{n_k_means}{(numeric) number of clusters to be identified when using the
"k-means" in \code{cluster_method}.}

\item{split_distance}{(numeric) distance in meters (if \code{space} = "G")
or Euclidean distance (if \code{space} = "E") to identify clusters if
\code{cluster_method} = "hierarchical".}
}
\value{
A data frame containing \code{data} and an additional column defining
clusters.
}
\description{
Finds clusters of data in two dimensions based on distinct
methods.
}
\details{
Clustering methods make distinct assumptions and one of them may perform
better than the other depending on the pattern of the data.

The k-means method tends to perform better when data are grouped spatially
(spherically) and clusters are of a similar size. The hierarchical
clustering algorithm usually takes more time than the k-means method. Both
methods make assumptions and may work well on some data sets but fail
on others.
}
\examples{
# Data
m_matrix <- read_master(system.file("extdata/m_matrix.rds",
                                    package = "biosurvey"))

# Cluster detection
clusters <-  find_clusters(m_matrix$data_matrix, x_column = "PC1",
                           y_column = "PC2", space = "E",
                           cluster_method = "hierarchical", n_k_means = NULL,
                           split_distance = 4)
head(clusters)
}
