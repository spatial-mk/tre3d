## Computation of crown complementarity (CC) between two tree point clouds
##
## Based on method of Williams et al., NEE, 2017, DOI: 10.1038/s41559-016-0063
##
## To assess the partitioning of canopy space, we introduce
## an index of crown complementarity that is based on the difference among trees
## in crown volume within strata from the ground to the top of the canopy.
## For each focal tree, crown volume was calculated from
## the crown dimensions measurements by treating each stratum
## as an elliptical cylinder.
## Crown complementarity (CC) was calculated for a pair of trees as
## the difference in crown volume (V) between the two trees (i and j) in each
## stratum (k) summed across all strata.
## The total difference in crown volume was expressed as a proportion
## of the combined volume of the pair of trees, as follows:
##
## CC = sum(abs(V_ik - V_jk)) / (V_i + V_j)
##
## Matthias Kunz, 26.02.2019
##
## Input: tree_i --> point cloud (x y z) of tree i
##        tree_j --> point cloud (x y z) of tree j

#' Computation of crown complementarity based on slice volumes.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Correction of QSM for multistem trees.
#' @references Williams et al., NEE, 2017, DOI: 10.1038/s41559-016-0063
#' @param tree_i A file or data.frame containing a tree point cloud with x y z
#' @param tree_j A file or data.frame containing a tree point cloud with x y z
#' @param strata_height Height of strata (in meter) for which the complementarity is computed. Default 0.3
#' @param vox Should the input point cloud be reduced to voxel grid for faster computation. Default to TRUE.
#' @param voxel_size Size (in meters) of the voxels when vox==TRUE. Default 0.03
#' @param plotting Should strata be plotted. Default to FALSE
#' @return The crown complementarity value (CC_ij)
#' @export
#' @examples
#' crown_complementarity(tree_i, tree_j, strata_height=0.3, plotting=FALSE)

## Definition of CC_ij function
crown_complementarity <- function(tree_i=NA, tree_j=NA , strata_height=0.3, vox=FALSE, voxel_size=0.03, plotting=FALSE){

  ## ----------------------------------
  ## Check what type of input is provided for tree_i: data.frame or file
  if(is.data.frame(tree_i)){
    # Leave as is
  } else{
    ## Check if file exists. If not, give warning
    if (!file.exists(tree_i)){
      stop(paste0("Point cloud file: ", tree_i, " does not exist."))
    } else {
      ## Read the file
      options(readr.num_columns = 0) ## suppress the output of readr
      tree_i <- readr::read_delim(file = tree_i, delim = " ", col_names = c("x","y","z"))
      }
  }
  ## Optional: Reduce to voxel grid
  if (vox==TRUE){
    tree_i <- cloud2vox(input_cloud = tree_i, voxel_size = voxel_size)
  }

  ## ----------------------------------
  ## Check what type of input is provided for tree_j: data.frame or file
  if(is.data.frame(tree_j)){
    # Leave as is
  } else{
    ## Check if file exists. If not, give warning
    if (!file.exists(tree_j)){
      stop(paste0("Point cloud file: ", tree_j, " does not exist."))
    } else {
      ## Read the file
      options(readr.num_columns = 0) ## suppress the output of readr
      tree_j <- readr::read_delim(file = tree_j, delim = " ", col_names = c("x","y","z"))
      }
  }
  ## Optional: Reduce to voxel grid
  if (vox==TRUE){
    tree_j <- cloud2vox(input_cloud = tree_j, voxel_size = voxel_size)
  }

  ## ----------------------------------
  ## Convert to data.frames
  tree_i <- as.data.frame(tree_i)
  tree_j <- as.data.frame(tree_j)

  ## ----------------------------------
  ## Get lowest point of the two trees
  z_min <- min(min(tree_i[,3]), min(tree_j[,3]))

  ## Get highest point of the two trees
  z_max <- max(max(tree_i[,3]), max(tree_j[,3]))

  ## Generate strata (30 cm) based on lowest point to highest point
  strata<-seq(from=z_min,to=z_max,by=strata_height)

  ## Add highest point strata if it is not in last strata
  if (z_max > max(strata)){strata=append(strata, max(strata)+strata_height)}

  ## Data frame that will store strata information
  strata_df <- data.frame(strata=strata, r_ik=0, r_jk=0, strata_height=strata_height)

  ## Set CC_ij, V_i, V_j, strata_sum
  CC_ij = 0
  V_i = 0
  V_j = 0
  strata_sum = 0

  ## Loop over strata and compute volume using convex hull
  for (k in 1:(length(strata)-1)){

    ## Select points in strata
    #cat("\nStrata",strata[k],"to", strata[k+1])
    strata_tree_i <- tree_i[which(tree_i[,3]<strata[k+1] & tree_i[,3]>strata[k]),]
    strata_tree_j <- tree_j[which(tree_j[,3]<strata[k+1] & tree_j[,3]>strata[k]),]

    ## Plotting
    if (plotting==T){
      if (nrow(strata_tree_i)>2){
        rgl::points3d(as.data.frame(strata_tree_i), col="grey", size=0.1)
        }
      if (nrow(strata_tree_j)>2){
        rgl::points3d(as.data.frame(strata_tree_j), col="grey", size=0.1)
        }
    }

    ## Convert to convex hull and then to polygon to extract area (if there are at least 3 points)
    ## Additionally, check if points are unique so they can form at least a triangle
    if (nrow(strata_tree_i)>2 & nrow(unique(strata_tree_i[,1:2]))>3){
      ## Compute convex hull
      strata_tree_i_ch <- chull(strata_tree_i[,1:2])
      coords_strata_i <- strata_tree_i[c(strata_tree_i_ch, strata_tree_i_ch[1]),1:2]  # closed polygon
      strata_tree_i_ch_area <- sp::Polygon(coords_strata_i)@area

      ## Compute theoretical cylinder radius in strata (r_ik)
      r_ik <- sqrt(strata_tree_i_ch_area/pi)
      strata_df[strata_df$strata==strata[k],]$r_ik = r_ik
      r = r_ik; n <- 300; theta <- seq(0, 2*pi, len=n)
      x <- cos(theta) * r
      y <- sin(theta) * r
      z <- rep(strata[k]+(strata_height/2), n)
      #lines3d(x,y,z, lwd=2, col="black")

      ## Plotting
      if (plotting==T){
        if (nrow(strata_tree_i)>2){
          coords_strata_j <- as.data.frame(coords_strata_i)
          rgl::lines3d(x=coords_strata_i[,1], y = coords_strata_i[,2], z=strata[k]+(strata_height/2), col="black", lwd=2)
          }
      }
    } else{
      strata_tree_i_ch_area = 0
    }

    ## Convert to convex hull and then to polygon to extract area (if there are at least 3 points)
    ## Additionally, check if points are unique so they can form at least a triangle
    if (nrow(strata_tree_j)>2 & nrow(unique(strata_tree_j[,1:2]))>3){
      ## Compute convex hull
      strata_tree_j_ch <- chull(strata_tree_j[,1:2])
      coords_strata_j <- strata_tree_j[c(strata_tree_j_ch, strata_tree_j_ch[1]),1:2]  # closed polygon
      strata_tree_j_ch_area <- sp::Polygon(coords_strata_j)@area

      ## Compute theoretical cylinder radius in strata (r_ik)
      r_jk <- sqrt(strata_tree_j_ch_area/pi)
      strata_df[strata_df$strata==strata[k],]$r_jk = r_jk
      r = r_jk; n <- 300; theta <- seq(0, 2*pi, len=n)
      x <- cos(theta) * r
      y <- sin(theta) * r
      z <- rep(strata[k]+(strata_height/2), n)
      #lines3d(x,y,z, lwd=2, col="blue")

      ## Plotting
      if (plotting==T){
        if (nrow(strata_tree_j)>2){
          coords_strata_j <- as.data.frame(coords_strata_j)
          rgl::lines3d(x=coords_strata_j[,1], y=coords_strata_j[,2], z=strata[k]+(strata_height/2), col="blue", lwd=2)
          }
      }
    } else{
      strata_tree_j_ch_area = 0
    }

    cat(nrow(strata_tree_i), " ", nrow(strata_tree_j),"\n")

    ## Convert to volume by multiplying with strata height
    V_ik = strata_tree_i_ch_area * strata_height
    V_jk = strata_tree_j_ch_area * strata_height
    strata_V_diff = abs(V_ik - V_jk)

    #cat("\nV_ik =",V_ik)
    #cat("\nV_jk =",V_jk,"\n")

    ## Sum the strata volumes and strata_diff
    strata_sum = strata_sum + strata_V_diff
    V_i = V_i + V_ik
    V_j = V_j + V_jk

  }

  ## Plot the data
  #x_min <- min(min(tree_i[,1]), min(tree_j[,1])); x_max <- max(max(tree_i[,1]), max(tree_j[,1])) # Get x extends
  #y_min <- min(min(tree_i[,2]), min(tree_j[,2])); y_max <- max(max(tree_i[,2]), max(tree_j[,2])) # Get y extends
  #plot(1,type="n", xlim=c(y_min, y_max), ylim=c(z_min, z_max), asp=1)
  #points(tree_i[,2:3], pch='.', col="black")
  #points(tree_j[,2:3], pch='.', col="blue")
  #for (k in 1:(length(strata))){abline(h = strata[k], col="red")}

  ## Generate plot for theoretical (cylindrical strata)
  if (plotting==T){
    strata_df$strata <- strata_df$strata - min(strata_df$strata)

    ## Compute GINI index as measure of inequality
    gini_i <- ineq::ineq(strata_df$r_ik,type="Gini")
    gini_j <- ineq::ineq(strata_df$r_jk,type="Gini")
    cat("\nGINI: i =", gini_i,"| j =", gini_j)

    strata_df_melt <- reshape2::melt(strata_df, id.vars=c("strata", "strata_height"))
    print(strata_df_melt)
    levels(strata_df_melt$variable) <- c("Tree i", "Tree j")
    strata_tree_i <- strata_df_melt[strata_df_melt$variable=="Tree i",]
    strata_tree_j <- strata_df_melt[strata_df_melt$variable=="Tree j",]
    p<-ggplot2::ggplot() +
      ggplot2::geom_bar(data=strata_tree_i, ggplot2::aes(x=strata+(strata_height/2), y=value, fill="red"),
               stat="identity", width = 0.3, color="black", alpha = 0.35) +
      ggplot2::geom_bar(data=strata_tree_j, ggplot2::aes(x=strata+(strata_height/2), y=value, fill="green"),
               stat="identity", width = 0.3, color="black", alpha = 0.35) +
      ggplot2::xlab("Height [m]") + ggplot2::ylab("Strata volume [m?]") +
      ggplot2::xlim(0, round(max(strata_df_melt$strata),0)+1) +
      ggplot2::ggtitle("Volume [m^3] per strata (strata height 30 cm)") +
      ggplot2::scale_fill_manual(name="", values=c("red","green"), breaks=c("red", "green"), labels=c("Tree i", "Tree j")) +
      ggplot2::scale_x_continuous(breaks = seq(0,20,2)) +
      ggplot2::theme_bw() +
      #facet_grid(.~variable) +
      ggplot2::annotate("text",0.3,max(strata_df_melt$value),
               label=paste0("CC =", round( (strata_sum / (V_i + V_j)),2)),
               hjust = 1.2) +
      ggplot2::coord_flip() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), panel.grid = ggplot2::element_blank(),
                     legend.position = c(0.9,0.9))
    print(p)
  }

  ## Compute the CC_ij and return
  CC_ij = strata_sum / (V_i + V_j)
  #cat("\nCC_ij", CC_ij)
  return(CC_ij)

} ## END-OF-FUNCTION
## ****************************************************************************
