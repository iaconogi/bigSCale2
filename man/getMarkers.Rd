% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellMethods.R
\name{getMarkers}
\alias{getMarkers}
\title{getMarkers}
\usage{
getMarkers(object)
}
\arguments{
\item{sce}{object of the SingleCellExperiment class.}
}
\value{
2D list, each row is a cluster, each column is a level.
}
\description{
Retrives a 2D list containing the markers of different specificity for each cluster.
For more information check the online tutorial at www.github.com
}
\examples{
Mlist=getMarkers(sce)

}
