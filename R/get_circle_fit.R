## ****************************************************************************
## Script that computes least squares fit of a circle into a set of 2d points (3d projected to xy)
## Matthias Kunz, 29.04.2019
## Note MK: try package 'conicfit'

#' Computes the 3d alpha shape, volume, surface area, and centroid of a 3D set of points
#' @author Matthias Kunz, last updated: 26.11.2019
#' @description computes the 2d alpha shape, 2d convex hull a 3D set of points as well as ground diameter and diameter at breast height\cr
#'     To avoid collinearity a small random jitter is added to the points before alpha-shape computation.
#' @param input_cloud A file or data.frame containing a tree point cloud with x y z. If file, it assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param plot_fit Should the circle fit be plotted. Default to FALSE.
#' @return A .txt files containing: CPA, CH, Height, GD, DBH, etc..
#' @examples
#' get_circle_fit(data.frame(x=runif(20),y=runif(20), z=runif(20)), plot=FALSE)
## ****************************************************************************

get_circle_fit <- function(input_cloud=data.frame(x=runif(20),y=runif(20),z=runif(20)), plot_fit=FALSE){

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
  ## Check if CloudCompare is installed
  if (file.exists("C:/Program Files/CloudCompare/CloudCompare.exe")==FALSE){
    stop(paste0("CloudCompare not installed under C:\\Progra~1\\CloudCompare\\CloudCompare.exe"))
  }

  current_wd <- getwd() # Remember original wd

  ## Write the data to a temp file
  input_temp <- paste0("C:/Temp/","temp_cloud_cirle_fit.xyz")
  write.table(df, input_temp, col.names = F, row.names = F, sep = " ")
  setwd(dirname(input_temp)) # Set wd

  ## Set parameter based on size of data
  npoints<- nrow(df)
  n_cc <- 40 #round(npoints/3,0)
  step_cc <- 5 # size for voxel grid reduction

  system("cmd.exe", input = paste0("C:\\Progra~1\\CloudCompare\\CloudCompare.exe", " -SILENT",
                                   " -O ", input_temp,
                                   " -AUTO_SAVE OFF", " -NO_TIMESTAMP",
                                   " -EXTRACT_CC ", step_cc ," ", n_cc,
                                   " -C_EXPORT_FMT",  " ASC", " -PREC 7", " -EXT", " xyz", " -SAVE_CLOUDS"))

  Sys.sleep(3) ## Wait until import is finished

  ## Helper function to make a data frame with points on a circle
  circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }

  cat("\nFit circle...\n")
  components<-list.files(path = "C:/Temp/", pattern = "COMPONENT*")
  if (length(components)==0){
    check_circle_fit<-data.frame(component=NA, r=NA, x=-9999, y=-9999, mean_dist=NA, size_check=NA)
  } else {
    check_circle_fit<-data.frame(component=components, r=NA, x=-9999, y=-9999, mean_dist=NA, size_check=NA)
  }

  ## Pruefe ob ueberhaupt componenten da sind
  if(length(components)==0){
    cat("\nOnly one component...\n")
    # wie vorher
    component<-read.table(input_temp, sep=" ")
    fit <- circular::lsfit.circle(component[,1], component[,2])
    r = fit$coefficients[1]
    x = fit$coefficients[2]
    y = fit$coefficients[3]

    check_circle_fit$component <- NA
    check_circle_fit$r <- r
    check_circle_fit$x <- x
    check_circle_fit$y <- y
    check_circle_fit$mean_dist <- NA
    check_circle_fit$size_check <- NA

    radius <- check_circle_fit

  } else {
    cat("\nMultiple components...\n")
    for (i in check_circle_fit$component){
      ## Fit circle for each component
      component <- read.table(i, sep = " ", header = )
      fit <- circular::lsfit.circle(component[,1], component[,2])
      r = fit$coefficients[1]
      x = fit$coefficients[2]
      y = fit$coefficients[3]

      check_circle_fit[check_circle_fit$component==i,]$r <- r
      check_circle_fit[check_circle_fit$component==i,]$x <- x
      check_circle_fit[check_circle_fit$component==i,]$y <- y

      ## Compute point to circle distance
      component$dist = NA
      for (ii in 1:nrow(component)){
        component[ii,]$dist <- abs(sqrt((component[ii,]$V1 - x)^2 + (component[ii,]$V2 -y )^2) - r)
      }
      mean_dist <- mean(component$dist, na.rm = T)
      check_circle_fit[check_circle_fit$component==i,]$mean_dist <- mean_dist

      ## Check if position and radius are within boundaries of GD points
      pxmax = max(component$V1)
      pymax = max(component$V2)
      pxmin = min(component$V1)
      pymin = min(component$V2)

      circle_points_radius = sqrt((pxmax - pxmin)^2 + (pymax - pymin)^2) / 2;
      if (circle_points_radius*1.2 < r){
        check_circle_fit[check_circle_fit$component==i,]$size_check <- "too large"
      } else {check_circle_fit[check_circle_fit$component==i,]$size_check <- "okay"}

      #  plot(component$V1, component$V2, asp=1, main=mean_dist)
      #  lines(circleFun(center = c(x,y), diameter = r*2))

    }

    cat("\nDelete temporary components")
    unlink(components)

    ## Pruefe ob nur eine connected component
    if (nrow(check_circle_fit)==1){
      radius <- check_circle_fit
    } else{
      check_circle_fit <- check_circle_fit[check_circle_fit$size_check=="okay",]
      radius<-check_circle_fit[check_circle_fit$mean_dist==min(check_circle_fit$mean_dist),]
    }
  }

  ## ----------------------------------
  ## Plot the fitted circle
  if(plot_fit==T){
    circle<-circleFun(center = c(radius$x, radius$y), diameter = radius$r*2)
    plot(df[,1], df[,2], pch='.', asp = 1)
    lines(circle, col="blue")
    points(radius$x, radius$y, col="blue", pch=7, cex=1.5)
    title(paste0("Radius [m]: ", round(radius$r,2)))
  }

  setwd(current_wd)
  cat("\n")
  print(radius)
  return(radius)
}
