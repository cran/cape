% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pheno_cor.R
\name{plot_pheno_cor}
\alias{plot_pheno_cor}
\title{Plot trait pairs against each other}
\usage{
plot_pheno_cor(
  data_obj,
  pheno_which = NULL,
  color_by = NULL,
  group_labels = NULL,
  text_cex = 1,
  pheno_labels = NULL,
  pt_cex = 1
)
}
\arguments{
\item{data_obj}{a \code{\link{Cape}} object}

\item{pheno_which}{A vector of trait names to plot. The default is to plot all traits.}

\item{color_by}{A character string indicating a value to color the traits by, for example sex or treatment.
It must be one of the covariates. See \code{\link{pheno2covar}}.}

\item{group_labels}{A vector of names for the legend indicating the groups for the colored dots.}

\item{text_cex}{A numeric value indicating the size of the text}

\item{pheno_labels}{A vector of names for traits to appear in the plot in case the column names are not very pretty.}

\item{pt_cex}{A numeric value indicating the size of the points.}
}
\description{
This function plots pairs of traits against each other
to visualize the correlations between traits.
}
\examples{
\dontrun{
plot_pheno_cor(data_obj, color_by = "sex", group_labels = c("Female", "Male"))
}

}
