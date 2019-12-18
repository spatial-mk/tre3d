## A function that computes the effective number of layers (ENL)
## See: Ehbrecht et al., 2017, Agricultural and Forest Meteorology, DOI:10.1016/j.agrformet.2017.04.012
## Matthias Kunz, 15.10.2019

#' Computes the effective number of layers (ENL) after Ehbrecht et al. 2016
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Computes the effective number of layers (ENL) after Ehbrecht et al. 2017
#' @param input_cloud A 3D point cloud file, e.g. "tree1.xyz".
#'     Assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param voxel_size Voxel size (in meter) for ENL computation. Default 0.1 m.
#' @param slice_thickness Thickness (in meter) of vertical layers for ENL computation. Default 0.5 m.
#' @return ENL value
#' @export
#' @examples
#' ENL(input_cloud=df, voxel_size=0.05, slice_thickness=0.5)

#library(VoxR) # 3D voxel grid

ENL <- function(input_cloud=df,
                 voxel_size=0.05, slice_thickness=0.5)
                 {

  cat("\nCompute Effectice Number of Layers (ENL)\n")
  ## Convert point cloud to voxel grid with side-length of voxel_size, default=0.05m
  voxel_size = voxel_size # default 0.05 meter
  df_vox <- VoxR::vox(input_cloud, res = voxel_size)

  ## Rename columns and use only first three columns
  df_vox <- df_vox[,1:3]
  colnames(df_vox) <- c("x","y","z")

  ## Set lowest tree point (along z axis) to zero
  df_vox$z <- df_vox$z - min(df_vox$z)

  ## Generate strata (default 0.50 m) based on lowest point to highest point
  slice_thickness = slice_thickness # default 0.5 meter
  strata <- seq(0, max(df_vox$z), slice_thickness)

  ## Generate data frame that holds result for each strata
  voxel_slices <- data.frame(lower=strata, upper=strata + slice_thickness,
                             voxel_in_slice=NA, proportion=NA)
  total_voxel = nrow(df_vox) ## Total number of voxels in all slices

  ## Loop over the strata
  for (ii in voxel_slices$lower){
    slice <- df_vox[df_vox$z >= ii & df_vox$z < ii + slice_thickness, ]
    voxel_in_strata <- nrow(slice)
    voxel_slices[voxel_slices$lower==ii,]$voxel_in_slice <- voxel_in_strata
    # Proportion of filled voxel in the strata to the sum of all voxels
    voxel_slices[voxel_slices$lower==ii,]$proportion <- voxel_in_strata / total_voxel
  }

  ## Compute ENL
  ENL = 1 / sum(voxel_slices$proportion^2, na.rm = T)
  result_ENL <- data.frame(ENL=ENL, voxel_size=voxel_size, slice_thickness=slice_thickness)

  rm(input_cloud)
  rm(df_vox)
  gc() # Garbage collection

  return(result_ENL)
}
