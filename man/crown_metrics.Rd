% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crown_metrics.R
\name{crown_metrics}
\alias{crown_metrics}
\title{Crown metrics.}
\usage{
crown_metrics(
  input_cloud = "*.xyz",
  slice_thickness = 0.3,
  cbh = 0,
  vox = TRUE,
  voxel_size = 0.03,
  use_concave_hull = FALSE,
  alpha_value = 1,
  plot_hulls = FALSE,
  pos_x = NA,
  pos_y = NA
)
}
\arguments{
\item{input_cloud}{A 3D point cloud file, e.g. "tree1.xyz".
Assumes first three columns to be x y z without header. Futher columns are ignored.}

\item{slice_thickness}{Thickness of vertical slices in meter. Default 0.3}

\item{cbh}{Crown base height in meter. Default 0.0}

\item{vox}{Should the input point cloud be reduced to voxel grid for faster computation. Default to TRUE.}

\item{voxel_size}{Size (in meters) of the voxels when vox==TRUE. Default 0.03}

\item{use_concave_hull}{Should slices be modelled as concave hulls instead of convex hull. Default to FALSE.}

\item{alpha_value}{Alpha value when use_convace_hull==TRUE. Default 1.0}

\item{plot_hulls}{Should the CPA and slice hull be plotted. Default to FALSE.}

\item{pos_x}{X coordinate of tree position. Default NA}

\item{pos_y}{Y coordinate of tree position. Default NA}
}
\value{
A data.frame with tree parameters
}
\description{
Computation of tree crown metrics.\cr
Currrently implemented: height, cbh, cpa cr_length,
    cr_length_to_heigh, cr_sinuosity, cr_compactnes, cr_density, cr_gini, cr_volume, cr_area cr_displacement
}
\examples{
crown_metrics("tree1.xyz", 0.25)
crown_metrics("tree1.xyz", 0.3, CBHfromQSM("tree1_qsm.txt"), 0.025)
}
\author{
Matthias Kunz, last updated: 27.04.2019
}
