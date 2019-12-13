## ****************************************************************************
## Script that computes the 3d alpha shape, volume, surface area, and centroid of a 3D set of points
## Matthias Kunz, 29.04.2019

#' Computes the 3d alpha shape, volume, surface area, and centroid of a 3D set of points
#' @author Matthias Kunz, last updated: 29.04.2019
#' @description computes the 3d alpha shape, volume, surface area, and centroid of a 3D set of points\cr
#'     To avoid collinearity a small random jitter is added to the points before alpha-shape computation.
#' @param input_cloud A file or data.frame containing a tree point cloud with x y z. If file, it assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param alpha Alpha value for alpha-shape computation. Default 0.5
#' @param plot Should the 3d alpha shape be plotted. Default to FALSE.
#' @return A data.frame containing: volume, surface area, vol_cen_x, vol_cen_y, vol_cen_z, alpha value and the alpha-shape object.
#' @export
#' @examples
#' get_3d_alpha_shape(data.frame(x=runif(20),y=runif(20), z=runif(20)), alpha=1.0)
## ****************************************************************************
get_3d_alpha_shape <- function(input_cloud=data.frame(x=runif(20),y=runif(20),z=runif(20)), alpha=1.0, plot=FALSE){

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
  ## Add random jitter to points
  df <- df + rnorm(nrow(df), mean=0.001, sd = 0.0005)

  ## ----------------------------------
  ## Compute the actual alphashape
  if (nrow(unique(df))>5){
    as3d <- alphashape3d::ashape3d(as.matrix(df[,1:3]), alpha=alpha, pert=T)
    volume <- alphashape3d::volume_ashape3d(as3d)
    kf <- alphashape3d::surfaceNormals(as3d)
    surface_area <- sum(sqrt(kf$normals[,1]^2 + kf$normals[,2]^2 + kf$normals[,3]^2))
  } else {volume=0; surface_area=0}

  ## Plotting
  if (plot==T){
    rgl::rgl.open()# Open a new RGL device
    rgl::rgl.bg(color = "white") # Setup the background color
    rgl::par3d(windowRect=c(50,50,1200,800)) # Change size of rgl window
    plot(as3d, col="darkgreen", transparency = 0.5,  triangles = TRUE, edges = FALSE, vertices = FALSE)
    #points3d(df, col="black")
    aspect3d("iso")
  }

  ## Return result
  result <- list("volume"=volume, "surface_area"=surface_area,
                 "vol_cen_x"=NA,  "vol_cen_y"=NA, "vol_cen_z"=NA, "alpha"=alpha,
                 "as3d"=as3d)
  return(result)

} ## END-OF-FUNCTION
## ****************************************************************************
