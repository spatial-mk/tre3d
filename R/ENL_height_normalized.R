## A function that computes the effective number of layers (ENL) for height normalized layers
## See: Ehbrecht et al., 2017, Agricultural and Forest Meteorology, DOI:10.1016/j.agrformet.2017.04.012
## Matthias Kunz, 27.02.2020

#' Computes the effective number of layers (ENL) after Ehbrecht et al. 2016
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Computes the effective number of layers (ENL) after Ehbrecht et al. 2017
#' @param input_cloud A 3D point cloud file as full path, e.g. "C:/tree1.xyz".
#'     Assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param voxel_size Voxel size (in meter) for ENL computation. Default 0.1 m.
#' @param slice_thickness Thickness (in meter) of vertical layers for ENL computation. Default 0.5 m.
#' @return ENL value
#' @export
#' @examples
#' ENL(input_cloud="*.xyz", voxel_size=0.05, slice_thickness=0.5)

#library(VoxR) # 3D voxel grid

ENL_height_normalized <- function(input_cloud=".xyz",
                 voxel_size=0.1, slice_thickness=0.5)
                 {

  ## Check if LAStools exist
  if (!dir.exists("C:/LAStools/")){stop("Please install LAStools under C:/LAStools")}

  ## Check if file exists
  if (!file.exists(input_cloud)){stop("File does not exist.")}

  ## Voxelisation for faster and simplyfied measurement
  cat(crayon::bold(crayon::blue("\nVoxelization\n")))
  voxelgrid <- paste0(dirname(input_cloud),"/", tools::file_path_sans_ext(basename(input_cloud)), "_vox.laz")
  system("cmd.exe", input = paste("C:\\LAStools\\bin\\lasvoxel.exe",
                                  "-i",input_cloud,
                                  "-o",voxelgrid,
                                  "-step 0.03 -v"))

  ## Classify ground points
  cat(crayon::bold(crayon::blue("\nClassify ground points\n")))
  classification <- paste0(dirname(input_cloud),"/", tools::file_path_sans_ext(basename(input_cloud)), "_class.laz")
  system("cmd.exe", input = paste("C:\\LAStools\\bin\\lasground.exe",
                                  "-i",voxelgrid,
                                  "-o",classification,
                                  "-all_returns -not_airborne -v -step 0.5",
                                  "-bulge 0.05 -spike 0.3 -coarse")) ## Step ggf. anpassen. Siehe lasground readme

  ## Compute height above ground (normalize height above ground) and classify slices at specific heights
  cat(crayon::bold(crayon::blue("\nNormalize height above ground\n")))
  height_norm <- paste0(dirname(input_cloud),"/", tools::file_path_sans_ext(basename(input_cloud)), "_height_norm.laz")
  system("cmd.exe", input = paste("C:\\LAStools\\bin\\lasheight.exe",
                                  "-i",classification,
                                  "-o",height_norm,
                                  "-cores 4", "-replace_z"))

  ## Extract the ground points
  cat(crayon::bold(crayon::blue("\nExtract the ground points\n")))
  ground <- paste0(dirname(input_cloud),"/", tools::file_path_sans_ext(basename(input_cloud)), "_ground.xyz")
  system("cmd.exe", input = paste("C:\\LAStools\\bin\\las2las.exe ",
                                   "-i ", classification,
                                   "-o ", ground,
                                   "-keep_class 2 -oparse xyz"))

  ## Extract height normalized points (without ground) as voxel grid
  cat(crayon::bold(crayon::blue("\nExtract height normalized points (excluding ground)\n")))
  voxelgrid_ENL_norm <- paste0(dirname(input_cloud),"/", tools::file_path_sans_ext(basename(input_cloud)), "_vox_ENL_norm.xyz")
  system("cmd.exe", input = paste("C:\\LAStools\\bin\\lasvoxel.exe",
                                  "-i",height_norm,
                                  "-o",voxelgrid_ENL_norm,
                                  "-drop_z_below ", voxel_size*2,
                                  "-step ", voxel_size,
                                  "-oparse xyz",
                                  " -v"))

  ## Extract height normalized ground voxel grid
  cat(crayon::bold(crayon::blue("\nExtract height normalized points (excluding ground)\n")))
  voxelgrid_ENL_norm_ground <- paste0(dirname(input_cloud),"/", tools::file_path_sans_ext(basename(input_cloud)), "_vox_ground_ENL_norm.xyz")
  system("cmd.exe", input = paste("C:\\LAStools\\bin\\lasvoxel.exe",
                                  "-i",height_norm,
                                  "-o",voxelgrid_ENL_norm_ground,
                                  "-oparse xyz",
                                  "-drop_z_above ", voxel_size*2,
                                  "-step ", voxel_size," -v"))


  cat(crayon::bold(crayon::blue("\nCompute Effectice Number of Layers (ENL)\n")))
  ## Convert point cloud to voxel grid with side-length of voxel_size, default=0.05m
  df_vox <- read.table(voxelgrid_ENL_norm, sep=" ", col.names = c("x","y","z"))

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

  ## Remove temporary files
  unlink(classification)
  unlink(height_norm)
  unlink(voxelgrid)
  unlink(ground)
  unlink(voxelgrid_ENL_norm)
  unlink(voxelgrid_ENL_norm_ground)

  rm(input_cloud)
  rm(df_vox)
  gc() # Garbage collection

  return(result_ENL)
}
