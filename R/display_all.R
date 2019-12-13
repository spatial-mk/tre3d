## ****************************************************************************
## Script that display a lot
## Matthias Kunz, 29.04.2019
## setwd("D:/Dropbox/TUD_KUNZ/tre3d/Code/tre3d/")
## devtools::document(); devtools::install(); devtools::load_all()
## tre3d::display_all(cloud2vox("D:/Dropbox/TUD_KUNZ/tre3d/Data/tree5.xyz"), cbh=1.4, qsm="D:/Dropbox/TUD_KUNZ/tre3d/Data/results/cyl_data_tree5.txt")
#' Script that display a lot
#' @author Matthias Kunz, last updated: 29.04.2019
#' @description Script that display a lot
#' @param input_cloud A file or data.frame containing a tree point cloud with x y z. If file, it assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param qsm qsm A file containing a QSM, e.g. cyl_data_*.txt
#' @param cbh Crown base height in meter. Default 0.0
#' @param alpha Alpha value for alpha-shape computation. Default 1.5
#' @param display_qsm Should all cylinders of the QSM be plotted. This may take a long time. Default to FALSE.
#' @param display_crown_volume Should 3d crown volume (alpha-shape) be plotted. Default to TRUE.
#' @param suppress_points Should 3d point be supressed from plotting. Default to FALSE.
#' @param n_seg Number of segment used for drawing cylinders. Default to 10.
#' @param random_subsampling Should the input point cloud be subsampled randomly. Default to TRUE.
#' @param mult Do you want to plot multiple trees to one device. Default to FALSE.
#' @param subsample_percent Numbers of points (in percent, min=5) after random subsampling. Default to 20 percent of original size.
#' @param set_to_zero Should the trees all be set to a elevation of zero. Ignores actual elevation. Default to TRUE.
#' @return Just plots stuff.
#' @examples
#' display_all(input_cloud="*.xyz", qsm="cyl_data.txt", cbh=0, alpha=1.0)
## ****************************************************************************
display_all <- function(input_cloud="*.xyz", random_subsampling=TRUE, qsm="cyl_data*.txt", cbh=0, alpha=1.5,
                        display_qsm=FALSE, display_crown_volume=TRUE, display_cpa=TRUE, display_strata=TRUE, suppress_points=FALSE,
                        n_seg=10, subsample_percent=25, set_to_zero=TRUE, mult=F){

  ## ----------------------------------
  ## Check what type of input is provided for tree_i: data.frame or file
  cat("\n\nGet the tree points")
  if(is.data.frame(input_cloud)){
    tree = input_cloud
  } else{
    ## Check if file exists. If not, give warning
    if (!file.exists(input_cloud)){
      stop(paste0("Point cloud file: ", input_cloud, " does not exist."))
    } else {
      ## Read the file, assuming column order: x y z
      options(readr.num_columns = 0) ## suppress the output of readr
      tree <- readr::read_delim(file = input_cloud, delim = " ", col_names = c("x","y","z"))
    }
  }

  ## Get the crown points
  z_offset <- min(tree[,3])
  tree_crown <- dplyr::filter(tree, z>z_offset+cbh)

  ## Set tree mininum height to zero
  if (set_to_zero==TRUE){
    tree[,3] <- tree[,3] - z_offset
    tree_crown[,3] <- tree_crown[,3] - z_offset
  }


  ## Subsample point cloud
  if (random_subsampling==TRUE){
    if (subsample_percent<5){subsample_percent=5} # Ensure a mininum of 5%
    cat("\nSubsample point cloud")
    tree<-dplyr::sample_n(tree, size=nrow(tree)*(subsample_percent/100))
  }

  ## ----------------------------------
  ## Read in QSM and check if file exists. If not, give warning
  if (display_qsm==TRUE){
    cat("\nGet the QSM")
    if (!file.exists(qsm)){
      stop(paste0("QSM file: ", qsm, " does not exist."))
    } else {
      ## Read the file, assuming column order: x y z
      options(readr.num_columns = 0) ## suppress the output of readr
      qsm_df <- read.table(file = qsm, sep= ",", header = T)
    }
    ## Set qsm mininum height to zero
    if (set_to_zero==TRUE){
      qsm_df[,5] <- qsm_df[,5] - min(qsm_df[,5])
    }

    ## Rename columns for better understanding
    qsm_df <- dplyr::rename(qsm_df, r=radius, l=length,
                            x=start_1, y=start_2, z=start_3,
                            xdir=axis_1, ydir=axis_2, zdir=axis_3)
  }

  ## ----------------------------------
  ## Get the crown projection area
  if (display_cpa==TRUE){
    cat("\nCompute CPA")
    cpa <- get_2d_alpha_shape(input_cloud= tree, alpha = alpha)
  }

  ## ----------------------------------
  ## Get the 3d crown alpha-shape
  if (display_crown_volume==TRUE){
    cat("\nCompute 3D (crown) alpha-shape")
    #tree_crown[,3] + z_offset
    as3d <- get_3d_alpha_shape(tree_crown, alpha = alpha)
  }

  ## ----------------------------------
  ## Plot everything with rgl
  clear_rgl=TRUE
  if (mult==F){
    while (rgl::rgl.cur() > 0) { rgl::rgl.close() }
    rgl::rgl.open()# Open a new RGL device
  } else {clear_rgl=FALSE}
  rgl::par3d(windowRect=c(50,50,1200,800), mouseMode = "trackball") # Change size of rgl window
  rgl::rgl.bg(color = "white") # Setup the background color

  ## Plot the 3d alpha-shape
  if (display_crown_volume==TRUE){
    cat("\nPlot 3D (crown) alpha-shape")
    plot(as3d$as3d, col="darkgreen", transparency = 0.25,  triangles = TRUE, edges = FALSE, vertices = FALSE, clear=clear_rgl)
  } else {}

  ## Plot the points
  if (suppress_points==FALSE){
    cat("\nPlot the tree points")
    tree[,3] <- tree[,3]
    rgl::points3d(tree, col="gray50")
  }

  ## ----------------------------------
  ## Get height strata
  if (display_strata==TRUE){
    strata_height=0.3
    cat("\nCompute strata")
    ## Get lowest point of the two trees
    z_min <- min(tree[,3])
    ## Get highest point of the two trees
    z_max <- max(tree[,3])
    ## Generate strata (30 cm) based on lowest point to highest point
    strata<-seq(from=z_min,to=z_max,by=strata_height)
    ## Add highest point strata if it is not in last strata
    if (z_max > max(strata)){strata=append(strata, max(strata)+strata_height)}
    ## Loop over strata and compute volume using convex hull
    for (k in 1:(length(strata)-1)){
      ## Select points in strata
      strata_tree <- tree[which(tree[,3]<strata[k+1] & tree[,3]>strata[k]),]

      ## Convert to convex hull and then to polygon to extract area (if there are at least 3 points)
      ## Additionally, check if points are unique so they can form at least a triangle
      if (nrow(strata_tree)>2 & nrow(unique(strata_tree[,1:2]))>3){
        ## Compute convex hull
        strata_tree_ch <- chull(strata_tree[,1:2])
        coords_strata <- strata_tree[c(strata_tree_ch, strata_tree_ch[1]),1:2]  # closed polygon
        ## Plotting
        coords_strata$z <- strata[k]+(strata_height/2)
        if (display_strata==TRUE){
          rgl::lines3d(coords_strata, col="black", lwd=2, add=T)
        }
      } else{}
    }
  }

  ## Plot the CPA
  if (display_cpa==TRUE){
    cat("\nPlot the CPA")
    if (set_to_zero==TRUE){cpa$polygon$z <- 0} else {cpa$polygon$z <- z_offset}
    rgl::polygon3d(cpa$polygon, col=rgb(1,1,0,1), alpha=0.5)
    rgl::lines3d(cpa$polygon, col="gray30")

    ## Plot the crown displacement (roughly) using tree position from lowest 10 cm slice
    cat("\nPlot the CPA displacement")
    if (set_to_zero==TRUE){
      elev = 0
      pos_x = mean(tree[tree$z<0.1,]$x)
      pos_y = mean(tree[tree$z<0.1,]$y)
    } else {
      elev=z_offset
      pos_x = mean(tree[tree$z<0.1+z_offset,]$x)
      pos_y = mean(tree[tree$z<0.1+z_offset,]$y)
    }
    rgl::arrow3d(c(pos_x, pos_y, elev),
                   c(cpa$cpa_cen_x,cpa$cpa_cen_y, elev),
                   type = "flat", col = "black", width = 0.2, thickness = 0.2, add=T)
    rgl::spheres3d(pos_x, pos_y, elev, col="black", r=0.02, add=T)
  }

  ## Plot the QSM
  # if (display_qsm == TRUE){
  #   cat("\nPlot the QSM")
  #   for (i in 1:nrow(qsm_df)){
  #     cylinder3d(p_start = c(qsm_df[i,c("x")], qsm_df[i,c("y")], qsm_df[i,c("z")]), # Start point of the cylinder
  #                orient = c(qsm_df[i,c("xdir")], qsm_df[i,c("ydir")], qsm_df[i,c("zdir")]), # Directional vector of the cylinder axis
  #                radius = qsm_df[i,c("r")],#/1000,      # Radius of the cylinder
  #                length = qsm_df[i,c("l")],#/1000,         # Length of the cylinder
  #                order = qsm_df[i,c("BranchOrder")],          # Branch order level of the cylinder
  #                n_seg = n_seg,         # Number of segments that draw the cylinder
  #                axis = FALSE,           # Display the cylinder axis
  #                cylinders = TRUE,      # Display the cylinder
  #                circles = FALSE,        # Display circles at the end of the cylinder
  #                transparent = TRUE)
  #   }
  # }

  ## Plot the QSM,  open3d()
  if (display_qsm == TRUE){

    cat("\nPlot the QSM")

    library(foreach)
    library(doParallel)
    cl <- makeCluster(2)
    registerDoParallel(cl)

    ## Generate branch order colors
    col_code = data.frame(BranchOrder=seq(from=0, to=9),
                          color=c("chocolate4","red","blue","darkgreen","darkorchid",
                                  "deeppink1","black","black","black","black"))
    qsm_df <- merge(qsm_df, col_code, by="BranchOrder", all.x=T)
    qsm_df$n_seg <- n_seg

    ## Plot the cylinders
    for (i in 1:nrow(qsm_df)){
      ## Get the cylinder axis
      p_start <- c(qsm_df[i,c("x")], qsm_df[i,c("y")], qsm_df[i,c("z")])
      orient = c(qsm_df[i,c("xdir")], qsm_df[i,c("ydir")], qsm_df[i,c("zdir")])
      p_end <- vector_end(p_start=p_start, orient=orient, length=qsm_df[i,c("l")])
      cylinder_axis <- matrix(c(p_start, p_end), ncol=3, byrow=T)
      ## Set color based on branch order
      #color = as.character(col_code[col_code == qsm_df[i,c("BranchOrder")],]$color)
      ## Plot the cylinder
      rgl::plot3d(rgl::cylinder3d(cylinder_axis, radius=qsm_df[i,c("r")],
                                  twist = 0, sides = qsm_df[i,c("n_seg")]), col=qsm_df[i,c("color")], alpha=0.5, add=T)
    }
  }

} ## END-OF-FUNCTION

#display_all(input_cloud = cloud2vox("D:/Dropbox/TUD_KUNZ/tre3d/Data/tree3.xyz", voxel_size=0.05),
#            qsm = "D:/Dropbox/TUD_KUNZ/tre3d/Data/results/cyl_data_tree3.txt",
#            alpha = 0.7,
#            cbh = cbh_from_qsm("D:/Dropbox/TUD_KUNZ/tre3d/Data/results/cyl_data_tree3.txt",
#                               min_diameter_m = 0.01, max_angle_deg = 80))

