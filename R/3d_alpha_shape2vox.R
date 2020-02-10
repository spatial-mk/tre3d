#' Voxelization of an 3D alphashape object
#' @author Carsten Hess, last updated: 10.02.2020
#' @description Voxelization of an 3D alpha shape crown model into the voxelgrid of specific voxel size
#' @param input_ashape3d Tree crowns ashape3d object based on alphashape3d R package (see get_3d_alpha_shape)
#' @param voxel_size Size of the voxels (in meter) of the underlying voxel grid. Default 0.03
#' @return Simple voxel grid point cloud (X, Y, Z) of all voxels inside the given 3D alpha shape.
#' @export
#' @examples
#'

ashape3D2vox <- function(input_ashape3d="", voxel_size=0.3) {

  # first simple check that an ashape3d object was delivered
  if(length(input_ashape3d)>0) {
    as3d <- input_ashape3d
  } else {
    stop('Missing ashape3d object! \n\n')
  }

  # get underlying point cloud of the ahsape3d
  pcd <- as.data.frame(as3d$x)

  # create closed bounding box voxelraster including the 3d alpha shape
  # estimate range dimensions of this bounding box
  xmin <- round_any(min(pcd$X), voxel_size, f=floor) - (2*voxel_size)
  xmax <- round_any(max(pcd$X), voxel_size, f=ceiling) + (2*voxel_size)
  xlen <- length(seq(xmin, xmax, voxel_size))

  ymin <- round_any(min(pcd$Y), voxel_size, f=floor) - (2*voxel_size)
  ymax <- round_any(max(pcd$Y), voxel_size, f=ceiling) + (2*voxel_size)
  ylen <- length(seq(ymin, ymax, voxel_size))

  zmin <- round_any(min(pcd$Z), voxel_size, f=floor) - (2*voxel_size)
  zmax <- round_any(max(pcd$Z), voxel_size, f=ceiling) + (2*voxel_size)
  zlen <- length(seq(zmin, zmax, voxel_size))

  # create data frame bounding box voxelraster as data frame
  x <- rep(seq(xmin, xmax, voxel_size), ylen*zlen)
  y <- rep(rep(seq(ymin, ymax, voxel_size), each=xlen), zlen)
  z <- rep(seq(zmin, zmax, voxel_size),each=xlen*ylen)

  # voxelraster of bounding box
  bbox <- data.frame(X=x, Y=y, Z=z)

  bbox.mat <- as.matrix(bbox)

  tryCatch({
    in3d <- inashape3d(as3d, indexAlpha=1, points = bbox.mat)
  }, warning = function(war) {
    in3d <- 0
  }, error = function(err) {
    in3d <- -1
  }, finally = {
    in3d <- inashape3d(as3d, indexAlpha=1, points = bbox.mat)
  })

  # Check IN/OUT TRUE/FALSE blue/red
  inoutcolors <- ifelse(in3d, "blue", "red")

  # Add inoutcolors to points
  bbox$cols <- inoutcolors

  # Subset of IN (blue) points
  bbox_in <- subset(bbox, cols=="blue", select=c(X,Y,Z))

  # fix possible rounding issues => to many digits
  bbox_in$X <- round_any(bbox_in$X, voxel_size)
  bbox_in$Y <- round_any(bbox_in$Y, voxel_size)
  bbox_in$Z <- round_any(bbox_in$Z, voxel_size)

  ## Voxelraster of a-shaped crown
  as3d_voxel <- unique(bbox_in)

  return(as3d_voxel)

}

