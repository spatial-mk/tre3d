% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSCI.R
\name{SSCI}
\alias{SSCI}
\title{SSCI computation after Ehbrecht et al. 2017}
\usage{
SSCI(
  input_fls = "*.fls",
  align_to_terrain = T,
  voxel_size = 0.1,
  slice_thickness = 0.5,
  angular_resolution = 360/2560,
  CROP_SOR_xyz = F,
  crop_square_length = 10,
  crop_circular = F,
  crop_circular_radius = 10,
  save_steps = F,
  plotting = F
)
}
\arguments{
\item{input_fls}{A .fls, .ptx or .xyz file. Note .fls requires a Main file in the same folder as the .fls file.}

\item{align_to_terrain}{Logical. Should a slope correction for SSCI be performed. Default TRUE.}

\item{voxel_size}{Voxel size (in meter) for ENL computation. Default 0.1 m.}

\item{slice_thickness}{Thickness (in meter) of vertical layers for ENL computation. Default 0.5 m.}

\item{angular_resolution}{Angular resolution of the scan, i.e. 360/number of points. Default 360/2560.}

\item{crop_square_length}{Half edge length (in meter) of squared crop along xy axis around scan center. Default +-10.}

\item{crop_circular}{Logical. Should points be clipped circular around the (0,0,0) position. Default FALSE.}

\item{crop_circular_radius}{Clip radius (in meter) if crop_circular=TRUE}

\item{save_steps}{Logical. Should aligned point clouds and filter clouds be saved. Default TRUE.}

\item{plotting}{Logical. Should the MeanFrac polygons be plotted. Note: This may takes a long time. Default FALSE.}

\item{CROP_SOR_XYZ}{Logical. Should a statistical outlier removal filter be applied and point cloud be cropped. Default FALSE.}
}
\value{
List with SSCI metrics.
}
\description{
Compute SSCI for .fls, .ptx or .xyz files
}
\examples{
SSCI()
}
\author{
Matthias Kunz, last updated: 18.12.2019
}
