% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cloud2vox.R
\name{cloud2vox}
\alias{cloud2vox}
\title{Voxel grid reduction of a point cloud.}
\usage{
cloud2vox(input_cloud = "*.xyz", voxel_size = 0.03, set_z_to_zero = FALSE)
}
\arguments{
\item{input_cloud}{A file or data.frame containing a tree point cloud with x y z. If file, it assumes first three columns to be x y z without header. Futher columns are ignored.}

\item{voxel_size}{Size of the voxels (in meter) in the reduced (gridded) point cloud. Default 0.03}

\item{set_z_to_zero}{Set lowest tree point (z axis) to zero. Default to FALSE.}
}
\value{
The reduced point cloud voxel grid as data.frame with column names x y z.
}
\description{
Converts point cloud (x y z) to a voxel grid with defined voxel size.
}
\examples{
cloud2vox("tree1.xyz", 0.03) # 0.7853982
}
\author{
Matthias Kunz, last updated: 27.04.2019
}
