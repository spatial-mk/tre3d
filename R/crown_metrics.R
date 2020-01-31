## ****************************************************************************
## An R-function that computes the following crown metrics from point clouds
##
## Crown Sinuosity --> see 10.1016/j.foreco.2016.04.047, FEM, 2016
## Crown Compactness --> MSc Fanghaenel, Circle--> (A/(U^2))*4*pi=1
## GINI index of crown volume --> library(ineq::ineq)
##
## Notes:
## + Variables require an estimate or measurement of crown base height (cbh)
##   If not available, set CBH to zero.
## + Script will generate slices of defined height (slice_thickness)
## + Slice hull is computed from convex hull
## + Point clouds are expected to have x y z values with space as delimeter,
##   with no header
## + For faster computation point cloud is reduced to voxel raster
##   with size voxel_size (default=2cm)
## + Input variables are in meter
##
## Script returns a data frame with tree name and the computed variables
## All variables are in meter or sqm or qbm
##
## Script requires the folling R-packages: dplyr, VoxR, ineq
##
## Author: Matthias Kunz, 04.04.2019
## ****************************************************************************

#' Crown metrics.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Computation of tree crown metrics.\cr
#' Currrently implemented: height, cbh, cpa cr_length,
#'     cr_length_to_heigh, cr_sinuosity, cr_compactnes, cr_density, cr_gini, cr_volume, cr_area cr_displacement
#' @param input_cloud A 3D point cloud file, e.g. "tree1.xyz".
#'     Assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param slice_thickness Thickness of vertical slices in meter. Default 0.3
#' @param cbh Crown base height in meter. Default 0.0
#' @param vox Should the input point cloud be reduced to voxel grid for faster computation. Default to TRUE.
#' @param voxel_size Size (in meters) of the voxels when vox==TRUE. Default 0.03
#' @param use_concave_hull Should slices be modelled as concave hulls instead of convex hull. Default to FALSE.
#' @param alpha_value Alpha value when use_convace_hull==TRUE. Default 1.0
#' @param plot_hulls Should the CPA and slice hull be plotted. Default to FALSE.
#' @param pos_x X coordinate of tree position. Default NA
#' @param pos_y Y coordinate of tree position. Default NA
#' @return A data.frame with tree parameters
#' @export
#' @examples
#' crown_metrics("tree1.xyz", 0.25)
#' crown_metrics("tree1.xyz", 0.3, CBHfromQSM("tree1_qsm.txt"), 0.025)

## The actual function
crown_metrics <- function(input_cloud="*.xyz",
                          slice_thickness=0.3, cbh=0, vox=TRUE, voxel_size=0.03,
                          use_concave_hull=FALSE, alpha_value=1.0, plot_hulls=FALSE,
                          pos_x=NA, pos_y=NA){

  ## ----------------------------------
  ## Check if file exists. If not, give warning
  if (!file.exists(input_cloud)){
    stop(paste0("Point cloud file: ", input_cloud, " does not exist."))
  }

  ## ----------------------------------
  ## Read in point cloud, assuming column order: x y z
  cat(paste0("Processing ",basename(input_cloud),": 0% -- "))
  ## Read in point cloud, assuming column order: x y z
  if (vox==TRUE){
    tree <- cloud2vox(input_cloud = input_cloud, voxel_size = voxel_size)
  } else {
    #tree <- read.table(input_cloud, sep=" ", header=F)
    #colnames(tree) <- c("x","y","z")
    options(readr.num_columns = 0) ## suppress the output of readr
    tree<-readr::read_delim(file = input_cloud, delim = " ", col_names = c("x","y","z"))
  }

  ## Set lowest tree point (z axis) to zero
  tree[,3] <- tree[,3] - min(tree[,3])

  ## ----------------------------------
  ## Compute tree height, crown length, and ratio
  height <- max(tree[,3]) # Height is in meter

  ## Check if CBH is greater than tree height from point cloud
  ## If so, set CBH to half of tree height to allow computation
  if (cbh > height){ cbh = height / 2
  warning("CBH > Height. Set CBH to half of tree height.")}

  cr_length = height - cbh
  cr_length_to_height <- cr_length / height

  ## ----------------------------------
  ## Use given tree position coordinates or
  ## get (rough) estimate of tree position based on stem base (0-15 cm)
  cat("25% -- ")
  if (is.na(pos_x) & is.na(pos_y)){
    stembase <- dplyr::filter(tree, z < 0.15)
    pos_x <- mean(stembase$x) ## Get tree position x coordinate
    pos_y <- mean(stembase$y) ## Get tree position y coordinate
  } else {
    pos_x=pos_x
    pos_y=pos_y
  }

  ## ----------------------------------
  ## Clip to crown (points above CBH) and set lowest crown z to zero
  crown <- dplyr::filter(tree, z > cbh)
  crown$z <- crown$z - min(crown$z)

  ## ----------------------------------
  ## Compute crown projection area and crown displacement
  ## Note: Adds small random noise so as3d is computed
  cat("50% -- ")
  if (nrow(unique(tree[,1:2]))>2){
    if (use_concave_hull==TRUE){
      cpa <- get_2d_alpha_shape(unique(tree[,1:2]), alpha = alpha_value, plot = plot_hulls)
      cpa$polygon <- NULL
      cpa <- as.data.frame(cpa)
    } else {
      cpa <- get_2d_chull(unique(tree[,1:2]), plot = plot_hulls)
      cpa$polygon <- NULL
      cpa <- as.data.frame(cpa)
    }
  } else {
    cpa=data.frame(area=0, perimeter=0, cpa_cen_x=pos_x, cpa_cen_y=pos_y)
  }

  cr_displacement <- sqrt((pos_x-cpa$cpa_cen_x)^2 + (pos_y-cpa$cpa_cen_y)^2)
  cr_disp_to_height <- cr_displacement / height
  cr_width_to_cr_length <- (sqrt(cpa$area/pi)) / cr_length

  ## ----------------------------------
  ## Compute crown sinuosity and crown compactness
  cat("60% -- ")

  ## Generate strata (30 cm) based on lowest point to highest point
  strata <- seq(0, cr_length, slice_thickness) # default 0.3 m slice

  ## Generate data frame that holds result for each strata
  result_slice <- data.frame(lower=strata, upper=strata + slice_thickness,
                             area=NA, perimeter=NA, radius=NA, compactness=NA)
  sinuosity=0 # set sinuosity to zero
  compactness=0 # set compactness to zero
  for (ii in result_slice$lower){
    schicht_centroid <- crown[crown$z >= ii & crown$z < ii + slice_thickness, ]
    schicht_centroid_xyz <- schicht_centroid[!duplicated(schicht_centroid[,1:2]),]
    if (nrow(unique(schicht_centroid_xyz))>3 & nrow(unique(schicht_centroid_xyz[,1:2]))>3){
      if (use_concave_hull==TRUE){
        ## Alpha_Hull_Variante
        center <- get_2d_alpha_shape(unique(schicht_centroid_xyz[,1:2]), alpha = alpha_value, plot=plot_hulls)
        center$polygon <- NULL
        center <- as.data.frame(center)
      } else {
        ## Convex_Hull_Variante
        center <- get_2d_chull(unique(schicht_centroid_xyz[,1:2]), plot=plot_hulls)
        center$polygon <- NULL
        center <- as.data.frame(center)
      }

    } else {
      ## Set area and perimeter to zero if there are less than three points in slice
      center=data.frame(area=0, perimeter=0, cpa_cen_x=pos_x, cpa_cen_y=pos_y)
    }
    slice_dist <- sqrt((pos_x-center$cpa_cen_x)^2 + (pos_y-center$cpa_cen_y)^2) # crown offset in slice
    compactness_slice <- (center$area / (center$perimeter^2)) * 4 * pi # compactness of the slice
    sinuosity=sinuosity+slice_dist # add distances to sinuosity per slice
    compactness = compactness + compactness_slice # add compactness of each slice
    r_ii <- sqrt(center$area/pi) ## Compute theoretical cylinder radius in slice (r_ii)
    result_slice[result_slice$lower==ii,]$area <- center$area
    result_slice[result_slice$lower==ii,]$perimeter <- center$perimeter
    result_slice[result_slice$lower==ii,]$radius <- r_ii
    result_slice[result_slice$lower==ii,]$compactness <- compactness_slice
  }
  #print(result_slice)
  sinuosity = sinuosity / cr_length
  compactness = mean(result_slice$compactness, na.rm = T)
  gini <- ineq::ineq(result_slice[result_slice$radius!=0,]$radius, type="Gini") # Compute GINI index

  ## ----------------------------------
  ## Crown volume and crown area from alpha-shape
  ## Note: Adds small random noise so as3d is computed
  cat("75% -- ")
  crown_as <- crown + rnorm(nrow(crown), mean=0.001, sd = 0.0002)

  ## Center crown points to 0,0,0 in case coordinates are too large, so that alphashape3d works
  crown_as[,1] <- crown_as[,1] - min(crown_as[,1])
  crown_as[,2] <- crown_as[,2] - min(crown_as[,2])
  crown_as[,3] <- crown_as[,3] - min(crown_as[,3])

  ## Check if enough crown points
  if (nrow(unique(crown_as))>5){
    as3d <- alphashape3d::ashape3d(as.matrix(crown_as[,1:3]), alpha=alpha_value, pert=T)
    cr_volume <- alphashape3d::volume_ashape3d(as3d)
    kf <- alphashape3d::surfaceNormals(as3d)
    cr_area <- sum(sqrt(kf$normals[,1]^2 + kf$normals[,2]^2 + kf$normals[,3]^2))
  } else {cr_volume=0; cr_area=0}
  cr_vol_to_cr_len = cr_volume / cr_length
  cr_area_to_cr_volume = cr_area / cr_volume
  cr_density = (nrow(unique(crown)) * (voxel_size^3)) / cr_volume

  #  library(rgl)
  #  plot(as3d, transparency=0.2, col="green", edges=F, vertices=F)
  #  points3d(crown, col="black", add=T)
  #  aspect3d("iso")
  #  bg3d("white")

  ## ----------------------------------
  ## Create a data frame that will hold the results
  metrics <- data.frame(tree=tools::file_path_sans_ext(basename(input_cloud)),
                        height=height,
                        cbh=cbh,
                        cpa=cpa$area,
                        cr_length=cr_length,
                        cr_length_to_height=cr_length_to_height,
                        cr_sinuosity=sinuosity,
                        cr_compactness=compactness,
                        cr_density=cr_density,
                        cr_gini=gini,
                        cr_volume=cr_volume,
                        cr_area=cr_area,
                        cr_displacement=cr_displacement,
                        cr_disp_to_height=cr_disp_to_height,
                        cr_vol_to_cr_len=cr_vol_to_cr_len,
                        cr_width_to_cr_length=cr_width_to_cr_length,
                        cr_area_to_cr_volume=cr_area_to_cr_volume)

  ## ----------------------------------
  ## Return the resulting data frame with the crown metrics
  cat("100%. Done.\n")
  return(metrics)

} # END-OF-FUNCTION
