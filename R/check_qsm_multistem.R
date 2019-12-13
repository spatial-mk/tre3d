## ****************************************************************************
#' Correction of QSM for multistem trees.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Correction of QSM for multistem trees.
#' @param qsm A file containing a QSM, e.g. cyl_data_*.txt
#' @param bifurcation_threshold Height threshold (in meter) under which bifurcation of an additional trunk is allowed. Default 1.5
#' @param percentage Percentage difference in radius between potential new trunk and main trunk. Default 40.0
#' @param diameter_treshold Minimum diameter (in meter) threhold for potential additional trunks. Default 0.035
#' @param min_tree_height Minimum tree height for searching for additional trunks. Default 2.0
#' @param plot_orders Plot branching orders original and checked QSM. Default to FALSE
#' @param plot_parts QSM for each detected stem (column tree_part). Default to FALSE
#' @param replace_qsm Should input QSM file be replaced. If FALSE _multistem.txt is added to file name. Default to FALSE
#' @return The path to the corrected QSM file.
#' @export
#' @examples
#' check_qsm_multistem("cyl_data*.txt", 1.5, 40.0, 0.035, 2.0, FALSE, FALSE, FALSE)
## ****************************************************************************
check_qsm_multistem <- function(qsm="cyl_data*.txt", bifurcation_threshold=1.5,
                                percentage=40.0, diameter_treshold=0.035, min_tree_height=2.0,
                                plot_orders=F, plot_parts=F, replace_qsm=F){
  ## ----------------------------------------------------------------------------
  ## Read in QSM
  cat("\nRead QSM", qsm)
  df <- read.table(qsm, header = T, sep = ",") # Original QSM

  ## ----------------------------------------------------------------------------
  ## Settings
  radius_percentage = percentage # Percentage of radius difference, default 40.0
  height_bifurcation = bifurcation_threshold # Max height [m] from which a trunk can form, default 1.5
  diameter_trunk = diameter_treshold # Minimum diameter [m] of potential new trunks
  min_tree_height=min_tree_height # Minimum height of tree that are checked for multistem

  offset=(max(df$start_1)-min(df$start_1)) + 1.5 # Only relevant for plotting

  ## ----------------------------------------------------------------------------
  ## Get lowest point and set to 0
  lowest <- dplyr::filter(df, BranchOrder==0, PositionInBranch==1)
  z_min <- min(lowest$start_3)
  df$start_3 <- df$start_3 - z_min

  ## Get tree height estimate
  height = max(df$start_3)

  ## Make row numbers as ID and create new data.frame for result
  df$ID <- seq.int(nrow(df))
  df$tree_part = 1 # as default only one part
  dfm <- df # dfm will hold the new/changed QSM

  ## Get all first order branches (only the start cylinder)
  pot_br <- dplyr::filter(df, BranchOrder==1, PositionInBranch==1)

  ## Only consider those that are larger than 3 cm in diameter
  pot_br <- dplyr::filter(pot_br, radius > (diameter_trunk/2))

  ## Only use those that start/bifurcate below 130 cm
  pot_br <- dplyr::filter(pot_br, start_3 < height_bifurcation)

  ## ----------------------------------------------------------------------------
  ## Only correct if there are candidates
  if (nrow(pot_br)!=0 & height>min_tree_height){
    cat("\nPotential multistem detected.")

    ## Add radius of parent cylinder at start
    cat("\nCheck radii.\n")
    pot_br$parent_radius = NA
    for (i in 1:nrow(pot_br)){
      parent_id <- pot_br[i,]$parent
      parent_radius <- df[df$ID==parent_id,]$radius
      pot_br[i,]$parent_radius <- parent_radius
    }

    ## ------------------------------------------------------------------------
    ## Only select candidates where circumreference difference is within treshold (radius_percentage) of parent
    #pot_br$r_diff = (abs(pot_br$radius - pot_br$parent_radius) / ((pot_br$radius + pot_br$parent_radius)/2)) * 100
    pot_br$r_diff = mapply(diff_percent, pot_br$radius, pot_br$parent_radius)
    cat("Candidates are:\n")
    print(dplyr::select(pot_br, ID, start_3, radius, parent_radius, r_diff))
    pot_br <- dplyr::filter(pot_br, r_diff < radius_percentage)
    cat(nrow(pot_br),"additional trunks detected.")

    new_part=1 # Counter for tree parts
    for (trunk in unique(pot_br$branch)){

      ## Get all cylinders of new trunks
      new_branches <- dplyr::filter(df, branch==pot_br[pot_br$branch==trunk,]$branch)
      dfm[dfm$branch%in%new_branches$branch,]$BranchOrder <- dfm[dfm$branch%in%new_branches$branch,]$BranchOrder - 1 # 1st to trunk
      dfm[dfm$branch%in%new_branches$branch,]$tree_part <- dfm[dfm$branch%in%new_branches$branch,]$tree_part + new_part

      for (order in 1:max(df$BranchOrder)){
        ## Get extentions of new trunks that are new first order branches
        potential_branches <- dplyr::filter(df, parent%in%new_branches$ID, BranchOrder>order)

        ## Get all cylinders new next order branches
        new_branches <- dplyr::filter(df, branch%in%unique(potential_branches$branch))

        ## Apply the changes
        dfm[dfm$branch%in%new_branches$branch,]$BranchOrder <- dfm[dfm$branch%in%new_branches$branch,]$BranchOrder - 1
        dfm[dfm$branch%in%new_branches$branch,]$tree_part <- dfm[dfm$branch%in%new_branches$branch,]$tree_part + new_part
      }

      new_part = new_part+1
    }

    ## ------------------------------------------------------------------------
    while (rgl::rgl.cur() > 0) { rgl::rgl.close() } # closes all open rgl windows

    ## Plotting of branches per branch order level
    if (plot_orders==T){
      ## Before and after correction
      rgl::rgl.open()# Open a new RGL device
      rgl::rgl.bg(color = "white") # Setup the background color
      rgl::par3d(windowRect=c(50,50,1200,800)) # Change size of rgl window
      rgl::plot3d(df$start_1,df$start_2,df$start_3, col="gray")
      rgl::clear3d()
      rgl::points3d(df$start_1,df$start_2,df$start_3, col="gray", box=F)
      rgl::text3d(min(df$start_1), min(df$start_2), min(df$start_3), texts = "QSM start points")
      rgl::text3d(min(df$start_1)+offset, min(df$start_2), min(df$start_3), texts = "QSM order before")
      rgl::text3d(min(df$start_1)+offset*2, min(df$start_2), min(df$start_3), texts = "QSM order multistem")
      rgl::aspect3d("iso")
      rgl::view3d(0, -60, zoom=1, fov=0)
      for (i in unique(df$branch)){
        # Original
        temp <- df[df$branch==i,]
        rgl::lines3d(temp$start_1+offset, temp$start_2, temp$start_3, lwd=2, col=temp$BranchOrder+1, add=T)

        # Corrected
        temp2 <- dfm[dfm$branch==i,]
        rgl::lines3d(temp2$start_1+offset*2, temp2$start_2, temp2$start_3, lwd=2, col=temp2$BranchOrder+1, add=T)
        #text3d(temp2[temp2$PositionInBranch==1,]$start_1+offset*2,
        #       temp2[temp2$PositionInBranch==1,]$start_2,
        #       temp2[temp2$PositionInBranch==1,]$start_3,
        #       texts = temp2$branch, cex=0.5, adj = 1.25)
      }
      rgl::axis3d('z', pos=c( 0, 0, NA ), col = "darkgrey")
    }

    ## Plotting of branches per tree part
    if (plot_parts==T){
      ## Before and after correction
      rgl::rgl.open()# Open a new RGL device
      rgl::rgl.bg(color = "white") # Setup the background color
      rgl::par3d(windowRect=c(50,50,1200,800)) # Change size of rgl window
      rgl::plot3d(df$start_1,df$start_2,df$start_3, col="gray")
      rgl::clear3d()
      rgl::points3d(df$start_1,df$start_2,df$start_3, col="gray", box=F)
      rgl::text3d(min(df$start_1), min(df$start_2), min(df$start_3), texts = "QSM start points")
      rgl::text3d(min(df$start_1)+offset, min(df$start_2), min(df$start_3), texts = "QSM order before")
      rgl::text3d(min(df$start_1)+offset*2, min(df$start_2), min(df$start_3), texts = "QSM order multistem")
      rgl::aspect3d("iso")
      rgl::view3d(0, -60, zoom=1, fov=0)
      for (i in unique(df$branch)){
        # Original
        temp <- df[df$branch==i,]
        rgl::lines3d(temp$start_1+offset, temp$start_2, temp$start_3, lwd=2, col=temp$tree_part, add=T)

        # Corrected
        temp2 <- dfm[dfm$branch==i,]
        rgl::lines3d(temp2$start_1+offset*2, temp2$start_2, temp2$start_3, lwd=2, col=temp2$tree_part, add=T)
        #text3d(temp2[temp2$PositionInBranch==1,]$start_1+offset*2,
        #       temp2[temp2$PositionInBranch==1,]$start_2,
        #       temp2[temp2$PositionInBranch==1,]$start_3,
        #       texts = temp2$branch, cex=0.5, adj = 1.25)
      }
      rgl::axis3d('z', pos=c( 0, 0, NA ), col = "darkgrey")
    }

    ## ------------------------------------------------------------------------
    ## Write result
    cat("\nQSM corrected.\n")
    dfm$ID <- NULL
    dfm$start_3 = dfm$start_3 + z_min ## Add height offset again

    if (replace_qsm==T){
      out_qsm=qsm # Replace old qsm file with updated version
    } else {
      out_qsm <- paste0(tools::file_path_sans_ext(qsm),"_multistem.txt") ## add multistem to output name
    }

    write.table(dfm, out_qsm, sep=",", col.names = T, row.names = F, quote = F)
  }  ## END OF MULTISTEM CASE
  else {
    cat("\nNo multistem detected.\n")
    df$ID <- NULL
    df$start_3 = df$start_3 + z_min ## Add height offset again
    write.table(df, qsm, sep=",", col.names = T, row.names = F, quote = F)
  }
  ## End-of-correction

  ## Return file path and name of corrected QSM
  return(out_qsm)

  cat("\n--------------")
} ## END-OF-FUNCTION
## ****************************************************************************
