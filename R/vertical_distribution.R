#' Vertical_distribution.
#' @author Matthias Kunz, last updated: 30.01.2020
#' @description Computation of vertical (wood) distribution
#' @param input_cloud A 3D point cloud file, e.g. "tree1.xyz".
#'     Assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param strata_height Thickness of vertical slices in meter. Default 0.3
#' @param vox Should the input point cloud be reduced to voxel grid for faster computation. Default to TRUE.
#' @param voxel_size Size (in meters) of the voxels when vox==TRUE. Default 0.03
#' @return A data.frame with tree parameters
#' @export
#' @examples
#' vertical_distribution("tree1.xyz", 0.25)

## The actual function
vertical_distribution <- function(input_cloud=data.frame(x=runif(20),y=runif(20),z=runif(20)), strata_height=0.3, vox=TRUE, voxel_size=0.03){

  ## ----------------------------------
  ## Check what type of input is provided for tree_i: data.frame or file
  if(is.data.frame(input_cloud)){
    df = input_cloud
  } else{
    ## Check if file exists. If not, give warning
    if (!file.exists(input_cloud)){
      stop(paste0("Point cloud file: ", input_cloud, " does not exist."))
    } else {
      ## Read the file, assuming column order: x y z
      options(readr.num_columns = 0) ## suppress the output of readr
      df <- readr::read_delim(file = input_cloud, delim = " ", col_names = c("x","y","z"))
    }
  }

  ## ----------------------------------
  ## Optional: Reduce to voxel grid
  if (vox==TRUE){
    df <- cloud2vox(input_cloud = df, voxel_size = voxel_size)
  }

  ## ----------------------------------
  ## Convert to data.frames
  df <- as.data.frame(df)

  ## ----------------------------------
  ## Set lowest tree point (z axis) to zero
  df[,3] <- df[,3] - min(df[,3])

  ## ----------------------------------
  ## Get lowest point of the two trees
  z_min <- min(df[,3])

  ## Get highest point of the two trees
  z_max <- max(df[,3])

  ## Generate strata (30 cm) based on lowest point to highest point
  strata<-seq(from=z_min,to=z_max,by=strata_height)

  ## Add highest point strata if it is not in last strata
  if (z_max > max(strata)){strata=append(strata, max(strata)+strata_height)}

  ## Data frame that will store strata information
  strata_df <- data.frame(strata=strata, n_points=0, r_k=0, strata_height=strata_height)

  ## Loop over strata and compute volume using convex hull
  for (k in 1:(length(strata)-1)){

    ## Select points in strata
    #cat("\nStrata",strata[k],"to", strata[k+1])
    strata_tree <- df[which(df[,3]<strata[k+1] & df[,3]>strata[k]),]

    ## Store number of points in strata
    strata_df[strata_df$strata==strata[k],]$n_points <- nrow(strata_tree)

    ## Convert to convex hull and then to polygon to extract area (if there are at least 3 points)
    ## Additionally, check if points are unique so they can form at least a triangle
    if (nrow(strata_tree)>2 & nrow(unique(strata_tree[,1:2]))>3){
      ## Compute convex hull
      strata_tree_ch <- chull(strata_tree[,1:2])
      coords_strata <- strata_tree[c(strata_tree_ch, strata_tree_ch[1]),1:2]  # closed polygon
      strata_tree_ch_area <- sp::Polygon(coords_strata)@area

      ## Compute theoretical cylinder radius in strata (r_ik)
      r_k <- sqrt(strata_tree_ch_area/pi)
      strata_df[strata_df$strata==strata[k],]$r_k = r_k
      r = r_k; n <- 300; theta <- seq(0, 2*pi, len=n)
      x <- cos(theta) * r
      y <- sin(theta) * r
      z <- rep(strata[k]+(strata_height/2), n)
      #lines3d(x,y,z, lwd=2, col="black")

    } else{
      strata_tree_ch_area = 0
    }

    ## Convert to volume by multiplying with strata height
    V_k = strata_tree_ch_area * strata_height

  } # end of for loop over strata

  return(strata_df)
}
