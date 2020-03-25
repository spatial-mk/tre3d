## ****************************************************************************
## Script that computes the 2d alpha shape, convex hull of the crown points
## as well as fits a circle at the GD and DBH height
## Matthias Kunz, 29.04.2019

#' Computes the 3d alpha shape, volume, surface area, and centroid of a 3D set of points
#' @author Matthias Kunz, last updated: 26.11.2019
#' @description computes the 2d alpha shape, 2d convex hull a 3D set of points as well as ground diameter and diameter at breast height\cr
#'     To avoid collinearity a small random jitter is added to the points before alpha-shape computation.
#' @param input_cloud A file or data.frame containing a tree point cloud with x y z. If file, it assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param vox Should the input point cloud be reduced to voxel grid for faster computation. Default to TRUE.
#' @param voxel_size Size (in meters) of the voxels when vox==TRUE. Default 0.03
#' @param alpha_value Alpha value. Default 1.0
#' @param overwrite_output Should outout files be overwritten?. Logical. Default FALSE.
#' @param plot Should the 2d alpha shape be plotted. Default to FALSE.
#' @param plot_hulls Should the CPA and slice hull be plotted. Default to FALSE.
#' @return A .txt files containing: CPA, CH, Height, GD, DBH, etc.
#' @export
#' @examples
#' get_3d_alpha_shape(data.frame(x=runif(20),y=runif(20), z=runif(20)), alpha=1.0)
## ****************************************************************************
single_tree_analysis <- function(input_cloud=data.frame(x=runif(20),y=runif(20),z=runif(20)),vox=FALSE, voxel_size=0.03, alpha_value=1.0,
                                 plot=FALSE, plot_hulls=FALSE, overwrite_output=FALSE){

  ## ----------------------------------
  ## Check if file was already processed and output exists
  file_id <- tools::file_path_sans_ext(basename(input_cloud))
  if (overwrite_output==FALSE & file.exists(paste0(dirname(input_cloud),"/analysis/",file_id,"_analysis.txt")) == TRUE){
    cat("\nOutput already exists. If you want to overwrite use overwrite_output==TRUE")
  }
  else if (overwrite_output==TRUE | file.exists(paste0(dirname(input_cloud),"/analysis/",file_id,"_analysis.txt")) == FALSE) {
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
        cat(paste0(crayon::bold(crayon::blue("\nProcessing ",basename(input_cloud))),":\n0% -- "))

        if (vox==TRUE){
          df <- cloud2vox(input_cloud = input_cloud, voxel_size = voxel_size)
        } else {
          options(readr.num_columns = 0) ## suppress the output of readr
          df<-readr::read_delim(file = input_cloud, delim = " ", col_names = c("x","y","z"))
        }

      }
    }

    ## ----------------------------------
    ## Set lowest tree point (z axis) to zero
    pos_low<-min(df[,3])
    df[,3] <- df[,3] - min(df[,3])

    ## ----------------------------------
    ## Get approximate position from lowest 5 cm in case fitting doesn,t work
    low_slice <- df[df[,3]<0.05,]
    pos_x <- mean(unlist(low_slice[,1]))
    pos_y <- mean(unlist(low_slice[,2]))

    ## ----------------------------------
    ## Compute tree height, crown length, and ratio
    height <- max(df[,3]) # Height is in meter

    ## ----------------------------------
    ## Add random jitter to points
    df <- df + rnorm(nrow(df), mean=0.001, sd = 0.0005)

    ## ----------------------------------
    ## Compute 2d alpha shape of all points projected to the xy plane
    cat("50% -- ")
    if (nrow(unique(df[,1:2]))>2){
      ## Alpha-shape
      cpa_as <- get_2d_alpha_shape(unique(df[,1:2]), alpha = alpha_value, plot = plot_hulls)
      cpa_as_poly <- cpa_as$polygon
      cpa_as$polygon <- NULL
      cpa_as <- as.data.frame(cpa_as)
      ## Convex hull
      cpa_ch <- get_2d_chull(unique(df[,1:2]), plot = plot_hulls)
      cpa_ch_poly <- cpa_ch$polygon
      cpa_ch$polygon <- NULL
      cpa_ch <- as.data.frame(cpa_ch)
    } else {
      cpa_as=data.frame(area=0, perimeter=0, cpa_cen_x=NA, cpa_cen_y=NA)
      cpa_ch=data.frame(area=0, perimeter=0, cpa_cen_x=NA, cpa_cen_y=NA)
    }

    ## ----------------------------------
    ## Compute ground diamater and diameter at breast height
    cat(paste0(crayon::bold(crayon::blue("\nCompute ground diameter and diameter at breast height\n"))))
    dbh_slice <- df[df[,3]<1.315 & df[,3]>1.295, ]
    gd_slice <- df[df[,3]<0.065 & df[,3]>0.045, ]

    ## For fit at least 10 points should exist, if not set to NA
    if (nrow(dbh_slice)<10){
      dbh_fit <- data.frame(component=NA, r=NA, x=-9999, y=-9999, mean_dist=NA, size_check=NA)
    } else {
      cat(crayon::green("\nGet diameter at breast height (130 cm)\n"))
      dbh_fit <- get_circle_fit(input_cloud = dbh_slice, plot_fit = TRUE)
    }

    ## Check if DBH fitting doesn't return anything
    if (nrow(dbh_fit)==0){
      dbh_fit <- data.frame(component=NA, r=NA, x=mean(dbh_slice[,1]), y=mean(dbh_slice[,2]), mean_dist=NA, size_check=NA)
    }

    if (nrow(gd_slice)<10){
      gd_fit <- data.frame(component=NA, r=NA, x=pos_x, y=pos_y, mean_dist=NA, size_check=NA)
    } else {
      cat(crayon::green("\nGet ground diameter (5 cm above ground)\n"))
      gd_fit <- get_circle_fit(input_cloud = gd_slice, plot_fit = TRUE)
    }

    ## Check if GD fitting doesn't return anything
    if (nrow(gd_fit)==0){
      gd_fit <- data.frame(component=NA, r=NA, x=mean(gd_slice[,1]), y=mean(gd_slice[,2]), mean_dist=NA, size_check=NA)
    }

    ## Check if GD is too large, i.e. radius > 100 cm,replace with mean xy of slice
    min_x <- min(gd_slice[,1]); max_x <- max(gd_slice[,1])
    min_y <- min(gd_slice[,2]); max_y <- max(gd_slice[,2])
    d_minx_maxx <- abs(max_x - min_x)
    d_miny_maxy <- abs(max_y - min_y)

    if (!is.na(gd_fit$r) & (gd_fit$r > ((d_miny_maxy+d_miny_maxy)/2))==T){
      cat("\nGD likely too big. Use mean values of slice coordinates instead.")
      print(gd_fit)
      print((d_miny_maxy+d_miny_maxy)/2)
      gd_fit$x <- pos_x
      gd_fit$y <- pos_y
      gd_fit$r <- sqrt((max_x-min_x)^2 + (max_y-min_y)^2) / 2
    }

    ## Write the result to analysis folder
    dir.create(file.path(dirname(input_cloud), "analysis"), showWarnings = FALSE)

    ## _ANALYSIS.txt
    if(overwrite_output==TRUE | file.exists(paste0(dirname(input_cloud),"/analysis/",file_id,"_analysis.txt")) == FALSE){
      sink(paste0(dirname(input_cloud),"/analysis/",file_id,"_analysis.txt"))

      cat(paste0("Tree_pos ", round(gd_fit$x,5), " ",round(gd_fit$y,5), " ",round(pos_low,5)," ", file_id," ", file_id, " ", file_id))
      cat(paste0("\nGD_m ", round(gd_fit$r*2,5), " ", file_id))
      cat(paste0("\nDBH_m ", round(dbh_fit$r*2,5), " ",round(dbh_fit$x,5)," ",round(dbh_fit$y,5)," ", file_id))
      cat(paste0("\nHeight_m ", round(height,5), " ", file_id))
      cat(paste0("\nCPA_CH_m ", round(cpa_ch$area,5), " ", file_id))
      cat(paste0("\nCentroid_CPA_CH ", round(cpa_ch$cpa_cen_x,5), " ",round(cpa_ch$cpa_cen_y,5), " ", file_id," ", file_id))
      cat(paste0("\nCPA_AS_m ", round(cpa_as$area,5), " ", file_id))
      cat(paste0("\nCentroid_CPA_AS ", round(cpa_as$cpa_cen_x,5), " ",round(cpa_as$cpa_cen_y,5), " ", file_id," ", file_id))
      cat(paste0("\nConvex_Hull_Polygon ",nrow(cpa_ch_poly)," ",file_id," ",file_id))
      for(i in 1:nrow(cpa_ch_poly)){
        cat(paste0("\n",round(cpa_ch_poly[i,]$x,5), " ", round(cpa_ch_poly[i,]$y,5), " ", file_id," ", file_id))
      }
      sink()
    }

    ## _CPA_AS.txt
    if(overwrite_output==TRUE | file.exists(paste0(dirname(input_cloud),"/analysis/",file_id,"_cpa_as.txt")) == FALSE){
      sink(paste0(dirname(input_cloud),"/analysis/",file_id,"_cpa_as.txt"))
      cat(paste0("Alpha_Shape_Polygon ",nrow(cpa_as_poly)," ",file_id," ",file_id))
      for(i in 1:nrow(cpa_as_poly)){
        cat(paste0("\n",round(cpa_as_poly[i,]$x,5), " ", round(cpa_as_poly[i,]$y,5), " ", file_id," ", file_id))
      }
      sink()
    }

    ## _DBH_CLOUD.xyz
    if(overwrite_output==TRUE | file.exists(paste0(dirname(input_cloud),"/analysis/",file_id,"_dbh_cloud.xyz")) == FALSE){
      dbh_slice[,3] <- dbh_slice[,3] + pos_low ## add height offset again
      write.table(dbh_slice, paste0(dirname(input_cloud),"/analysis/",file_id,"_dbh_cloud.xyz"), row.names = F, col.names = F, sep = " ", quote = F)
    }

    ## _GD_CLOUD.xyz
    if(overwrite_output==TRUE | file.exists(paste0(dirname(input_cloud),"/analysis/",file_id,"_gd_cloud.xyz")) == FALSE){
      gd_slice[,3] <- gd_slice[,3] + pos_low ## add height offset again
      write.table(gd_slice, paste0(dirname(input_cloud),"/analysis/",file_id,"_gd_cloud.xyz"), row.names = F, col.names = F, sep = " ", quote = F)
    }

    ## Print the result as a PNG to check
    cat(crayon::bold("\nAnalysis done. Print PNG\n"))
    tre3d::print_CPA(input_cloud)
  }

} ## END-OF-FUNCTION
