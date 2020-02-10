#' Voxel grid reduction an interpolation of ground point cloud (DGM).
#' @author Carsten Hess, last updated: 10.02.2020
#' @description Converts cloud of ground points (x y z) to a voxel grid with defined voxel size.
#' @param input_dgm A file or data.frame containing a ground points with x y z. If file, it assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param voxel_size Size of the voxels (in meter) in the reduced (gridded) point cloud of the ground (DGM). Default 0.1
#' @param raster_size Size of the rastergrid (in meter) to preprocess the ground points (DGM) to further eliminate noise points above ground.
#' @return The reduced point cloud voxel grid as data.frame with column names x y z.
#' @export
#' @examples
#' cloud2vox("tree1.xyz", 0.03) # 0.7853982
get_dgm_voxelmodel <- function(input_cloud="*.xyz", voxel_size=0.1, raster_size=0.2){

  ## Check what type of input is provided for tree_i: data.frame or file
  if(is.data.frame(input_cloud)){
    pcd = input_cloud
  } else{
    ## Check if file exists. If not, give warning
    if (!file.exists(input_cloud)){
      stop(paste0("Point cloud file: ", input_cloud, " does not exist."))
    } else {
      ## Read the file, assuming column order: x y z
      options(readr.num_columns = 0) ## suppress the output of readr
      pcd <- readr::read_delim(file = input_cloud, delim = " ")
    }
  }

  dgm <- pcd[,c(1:3)]
  colnames(dgm)<-c("X","Y","Z")

  ## rasterize high resolution DGM to eliminate possible noise
  ## raster_size > voxel_size
  raster <- dgm
  raster$X <- round_any(raster$X, raster_size)
  raster$Y <- round_any(raster$Y, raster_size)

  ## dummy tables for data manipulation
  dy <- aggregate(Z ~ X + Y, raster, min)

  dy1 <- interp(dy$X, dy$Y, dy$Z, duplicate="mean"
                , xo=seq(min(dy$X), max(dy$X), by = voxel_size)
                , yo=seq(min(dy$Y), max(dy$Y), by = voxel_size)
  )

  ## convert interpolation into a data.frame
  dy2 <- data.frame(X=rep(dy1$x,length(dy1$y))
                    , Y=rep(dy1$y,each=length(dy1$x))
                    , Z=as.vector(dy1$z))

  ## remove NAs
  dy2 <- na.omit(dy2)

  ## transform Z into the voxelraster
  dy2$Z <- round_any(dy2$Z, voxel_size)

  ## fix roundings
  dy2$X <- round_any(dy2$X, voxel_size)
  dy2$Y <- round_any(dy2$Y, voxel_size)

  ## Return DGM of unique ground voxels
  return(unique(subset(dy2, select=c(X,Y,Z))))

} # END-OF-FUNCTION
