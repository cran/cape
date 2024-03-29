% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_full_network.R
\name{plot_full_network}
\alias{plot_full_network}
\title{Plot the final epistatic network in a traditional network view.}
\usage{
plot_full_network(
  data_obj,
  p_or_q = 0.05,
  collapsed_net = TRUE,
  main = NULL,
  color_scheme = c("DO/CC", "other"),
  pos_col = "brown",
  neg_col = "blue",
  bg_col = "gray",
  light_dark = "f",
  node_border_lwd = 1,
  layout_matrix = NULL,
  zoom = 1,
  xshift = 0,
  yshift = 0,
  node_radius = 1,
  label_nodes = TRUE,
  label_offset = 0,
  label_cex = 1,
  legend_radius = 1,
  legend_cex = 1,
  legend_position = "topleft",
  arrow_offset = node_radius,
  arrow_length = 0.2,
  edge_lwd = 2
)
}
\arguments{
\item{data_obj}{A \code{\link{Cape}} object}

\item{p_or_q}{The maximum p value (or q value if FDR was used) for significant 
main effects and interactions.}

\item{collapsed_net}{A logical value indicating whether to show the network
condensed into linkage blocks (TRUE) or each individual marker (FALSE)}

\item{main}{A title for the plot}

\item{color_scheme}{either "DO/CC" or "other". "DO/CC" uses the official "DO/CC"
colors for the DO/CC alleles  
\url{http://www.csbio.unc.edu/CCstatus/index.py?run=AvailableLines.information}
"other" uses an unrelated color palette for multiple alleles.}

\item{pos_col}{The color to use for positive main effects and interactions
must be one of "green", "purple", "red", "orange", "blue", "brown", "yellow", "gray"
see \code{\link{get_color}}}

\item{neg_col}{The color to use for negative main effects and interactions
takes the same values as pos_col.}

\item{bg_col}{The color to be used in pie charts for non-significant main effects.
Takes the same values as pos_col}

\item{light_dark}{Indicates whether pos_col, neg_col, and bg_col should be selected
from light colors ("l"), dark colors ("d") or the full spectrum from light to dark ("f")}

\item{node_border_lwd}{The thickness of the lines around the pie charts}

\item{layout_matrix}{Users have the option of providing their own layout matrix for the
network. This should be a two column matrix indicating the x and y coordinates of each 
node in the network.}

\item{zoom}{Allows the user to zoom in and out on the image if the network is either 
running off the edges of the plot or too small in the middle of the plot.}

\item{xshift}{A constant by which to shift the x values of all nodes in the network.}

\item{yshift}{A constant by which to shift the y values of all nodes in the network.}

\item{node_radius}{The size of the pie chart for each node.}

\item{label_nodes}{A logical value indicating whether the nodes should be labeled.
Users may want to remove labels for large networks.}

\item{label_offset}{The amount by which to offset the node labels from the center of
the nodes.}

\item{label_cex}{The size of the node labels}

\item{legend_radius}{The size of the legend indicating which pie piece corresponds to which
traits.}

\item{legend_cex}{The size of the labels in the legend.}

\item{legend_position}{The position of the legend on the plot}

\item{arrow_offset}{The distance from the center of the node to the arrow head.}

\item{arrow_length}{The length of the head of the arrow}

\item{edge_lwd}{The thickness of the arrows showing the interactions.}
}
\value{
This function invisibly returns a list of length two. The first element
contains the igraph network object. The second contains the layout matrix for the
network. This can be passed in as an argument ("layout_matrix") which provides more
control to the user in the layout. Other network layouts from igraph can also be passed 
in here.
}
\description{
This function plots the final results in a layout different to
 both \code{\link{plot_variant_influences}} and \code{\link{plot_network}}. 
 In this view, the network is plotted with a traditional network layout. 
 The genomic position information in \code{\link{plot_network}} is lost, but 
 in this view it is easier to see the structure of the overall network 
 in terms of hubs and peripheral nodes. In this view, each node is plotted 
 as a pie-chart, and the main effects of the node are indicated as 
 positive, negative, or not-significant (gray). Significant 
 interactions are shown arrows between 
 nodes and colored based on whether they are positive or negative interactions. 
 Colors for positive and negative main effects and interactions are specified
in the arguments. The function \code{\link{get_network}} must be run before plotting 
 the network.
}
\details{
For most networks, the default options will be fine, but there is a lot of room
for modification if changes are desired
}
\references{
Csardi G, Nepusz T: The igraph software package for complex network 
research, InterJournal, Complex Systems 1695. 2006. \url{https://igraph.org/}
}
