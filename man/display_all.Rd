% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display_all.R
\name{display_all}
\alias{display_all}
\title{Script that display a lot}
\usage{
display_all(
  input_cloud = "*.xyz",
  random_subsampling = TRUE,
  qsm = "cyl_data*.txt",
  cbh = 0,
  alpha = 1.5,
  display_qsm = FALSE,
  display_crown_volume = TRUE,
  display_cpa = TRUE,
  display_strata = TRUE,
  suppress_points = FALSE,
  n_seg = 10,
  subsample_percent = 25,
  set_to_zero = TRUE,
  mult = F
)
}
\arguments{
\item{input_cloud}{A file or data.frame containing a tree point cloud with x y z. If file, it assumes first three columns to be x y z without header. Futher columns are ignored.}

\item{random_subsampling}{Should the input point cloud be subsampled randomly. Default to TRUE.}

\item{qsm}{qsm A file containing a QSM, e.g. cyl_data_*.txt}

\item{cbh}{Crown base height in meter. Default 0.0}

\item{alpha}{Alpha value for alpha-shape computation. Default 1.5}

\item{display_qsm}{Should all cylinders of the QSM be plotted. This may take a long time. Default to FALSE.}

\item{display_crown_volume}{Should 3d crown volume (alpha-shape) be plotted. Default to TRUE.}

\item{suppress_points}{Should 3d point be supressed from plotting. Default to FALSE.}

\item{n_seg}{Number of segment used for drawing cylinders. Default to 10.}

\item{subsample_percent}{Numbers of points (in percent, min=5) after random subsampling. Default to 20 percent of original size.}

\item{set_to_zero}{Should the trees all be set to a elevation of zero. Ignores actual elevation. Default to TRUE.}

\item{mult}{Do you want to plot multiple trees to one device. Default to FALSE.}
}
\value{
Just plots stuff.
}
\description{
Script that display a lot
}
\examples{
display_all(input_cloud="*.xyz", qsm="cyl_data.txt", cbh=0, alpha=1.0)
}
\author{
Matthias Kunz, last updated: 29.04.2019
}
