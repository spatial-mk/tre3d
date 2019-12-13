#' Display tree crown base height.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Display tree crown base height.
#' @param input_cloud A 3D point cloud file, e.g. "tree1.xyz".
#'     Assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param cbh Crown base height in meter. Default 0.0
#' @param vox Should the input point cloud be displayed as voxel grid for faster display. Default to TRUE.
#' @param voxel_size Size (in meters) of the voxels when vox==TRUE. Default 0.03
#' @return NULL. Function currently does not return anything.
#' @export
#' @examples
#' display_cbh("tree1.xyz", 3.4)
display_cbh <- function(input_cloud="*.xyz", cbh=0, vox=TRUE, voxel_size=0.03){

  ## Check if file exists. If not, give warning
  if (!file.exists(input_cloud)){
    stop(paste0("Point cloud file: ", input_cloud, " does not exist."))
  }

  ## Read in point cloud, assuming column order: x y z
  if (vox==TRUE){
    tree <- cloud2vox(input_cloud)
  } else {
    tree <- read.table(input_cloud, sep=" ", header=F)
  }

  ## Set lowest tree point (z axis) to zero
  tree[,3] <- tree[,3] - min(tree[,3])

  ## Plotting
  while (rgl::rgl.cur() > 0) { rgl::rgl.close() } # closes all open rgl windows
  rgl::rgl.open()# Open a new RGL device
  rgl::rgl.bg(color = "white") # Setup the background color
  rgl::plot3d(tree, col="gray10")
  rgl::planes3d(a=0,b=0,c=-1,d=cbh, col="red", alpha=0.5)
  rgl::aspect3d("iso")
  rgl::par3d(windowRect=c(50,50,1200,800), mouseMode = "trackball") # Change size of rgl window

} ## END-OF-FUNCTION
