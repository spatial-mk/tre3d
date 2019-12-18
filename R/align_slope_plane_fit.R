## ****************************************************************************
## A function that takes a 3d point cloud fits a plane using SVD and rotates
## Matthias Kunz, 15.08.2019
## Input: Terrain data
## Output: Rotation matrix
## Note: assummes 0,0,0 as rotation center, so reduced coordinates to the mean
##       then apply rotation R via: as.matrix(df[,1:3]) %*% R
##       then add the mean coordinates again
## ****************************************************************************

#' 3D plane fit for a set of points
#' @author Matthias Kunz, last updated: 18.12.2019
#' @description Fits a 3D plane into a set of 3D points (x,y,z)
#' @param terrain_data A file containing the (terrain) points into which the plane will be fitted
#' @param full_data A file containing the points that will be aligned/rotation based on the plane surface normal
#' @return Deos not return anything at the moment but produces files.
#' @export
#' @examples
#' align_slope_plane_fit(terrain_data="terrain.xyz", full_data="full.xyz", plotting=F)

## Load libraries
#library(readr) # Fast file input
#library(expm) # Compute cross-product

## ++++++
## INPUTS
#terrain <- read.table("D:/Kunz/Desktop/SSCI_Test/SSCI_E34.fls/SSCI_E34_ground_dem.xyz", sep=" ", col.names = c("x","y","z")) #c("x","y","z","x1","y1","z1","i")
#full_dataset <- read_delim("D:/Kunz/Desktop/SSCI_Test/SSCI_E34.fls/SSCI_E34_CROPPED_SOR.xyz", delim = " ",col_names = F)
#output <- "D:/Kunz/Desktop/SSCI_Test/SSCI_E34.fls/SSCI_E34_CROPPED_SOR_aligned.xyz"
#rotation_matrix <- "D:/Kunz/Desktop/SSCI_Test/SSCI_E34.fls/cc_matrix.txt"
## ++++++

align_slope_plane_fit <- function(terrain_data="*.xyz", full_data="*.xyz", plotting=F){
  current_wd <- getwd() # Remember original wd
  setwd(dirname(full_data)) # Set wd

  ## ----------------------------------------------------------------------------
  ## Step 0: Read in the data
  terrain <- suppressMessages(readr::read_delim(terrain_data, delim=" ", col_names = F))
  full_dataset <- suppressMessages(readr::read_delim(full_data, delim = " ",col_names = F))
  rotation_matrix <- paste0(tools::file_path_sans_ext(basename(full_data)),"_cc_matrix.txt")
  output <- paste0(tools::file_path_sans_ext(basename(full_data)),"_aligned.xyz")
  output2 <- paste0(tools::file_path_sans_ext(basename(terrain_data)),"_aligned.xyz")

  ## ----------------------------------------------------------------------------
  ## Step 1: Read in (terrain) points and reduced to mean coordinates
  #df <- read.table("D:/Darien_12_las/results/plot_cloud_lowest.xyz", sep=",", col.names = c("x","y","z"))
  #df <- read.table("W:/Alexandra/Darien_plot_04_las/results/plot_cloud_lowest.xyz", sep=",", col.names = c("x","y","z"))
  ##df <- read.table("D:/Kunz/Desktop/ground_points.xyz", sep=" ", col.names = c("x","y","z","x1","y1","z1","i"))
  ## Reduce coordinates to mean of the data set
  df <- as.data.frame(terrain)
  colnames(df) <- c("x","y","z")
  mean_x <- mean(df$x)
  mean_y <- mean(df$y)
  mean_z <- mean(df$z)
  reduced <- df[,1:3]
  reduced$x <- reduced$x - mean_x
  reduced$y <- reduced$y - mean_y
  reduced$z <- reduced$z - mean_z

  cat("X offset:",mean_x," | Y offset:",mean_y," | Z offset:",mean_z)

  ## ----------------------------------------------------------------------------
  ## Step: Fit a plane using SVD
  cat(bold(blue("\n>>> Compute plane fit using SVD")))
  ## https://math.stackexchange.com/questions/2810048/plane-fitting-using-svd-normal-vector
  ## Compute SVD and take v3 from V matrix to get plane normal vector
  SVD <- svd(reduced)
  a <- SVD$v[1,3]
  b <- SVD$v[2,3]
  c <- SVD$v[3,3]
  ## Check if surface plane vector is pointing the same direction (upwards) as the z axis
  if (c<0){a=-a;b=-b;c=-c}
  ## Compute d
  d <- -(a*mean_x + b*mean_y + c*mean_z) # get d of plane equation using one point on the plane

  ## Compute signed distance of all points to the plane
  ## Function to compute signed plane distance
  ## http://mathworld.wolfram.com/Point-PlaneDistance.html
  signed_dist2plane <- function(p=c(0,0,0), a=1, b=1, c=0, d=0){
    x0 = p[1]; y0 = p[2]; z0 = p[3]
    D = (a*x0 + b*y0 + c*z0 + d) / sqrt(a^2 + b^2 + c^2)
    if(is.na(D)){D=0}
    return(D)
  }
  df$dist <- NA
  for (i in 1:nrow(df)){
    p <- c(df[i,1], df[i,2], df[i,3])
    df[i,]$dist <- (signed_dist2plane(p=p, a=a, b=b, c=c, d=d))
  }
  ## Print and plot the points and the fitted plane
  cat("\n\na=", a, " b=", b, " c= ", c, " d=", d, " (meanErr=",mean(df$dist, na.rm=T),")", sep="")

  ## ----------------------------------------------------------------------------
  ## Step 3: Rotate coordinates to make ground paralell to plane (ax + by + cz + d = 0)
  cat(crayon::bold(crayon::blue("\n>>> Rotate coordinates\n")))
  ## https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
  rotate2vertical <- function(direction){

    ## Make normal vector to unit vector
    direction <- direction / sqrt(direction[1]^2 + direction[2]^2 + direction[3]^2)

    ## Unit vector along the z-axis
    vertical_vector <- c(0,0,1)

    ## Define vector for direction (b) and vertical vector (a)
    ## Note that a and b are exchanged compared to the solution of stackexchange
    a <- vertical_vector
    b <- direction

    ## Cross product of a and b
    v <- c(a[2]*b[3] - a[3]*b[2],
           a[3]*b[1] - a[1]*b[3],
           a[1]*b[2] - a[2]*b[1])

    ## Identity matrix
    I <- matrix(c(1,0,0,
                  0,1,0,
                  0,0,1), nrow=3, byrow=T)

    # Magitude/Length of v
    s <- sqrt(v[1]^2 + v[2]^2 + v[3]^2)

    ## Dot product of a and b
    c = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

    ## Skew-symetric cross-product of v
    v_cross <- matrix(c(0,-v[3],v[2],
                        v[3],0,-v[1],
                        -v[2],v[1],0), nrow=3, byrow=T)

    ## Compute rotation matrix (library(expm))
    R = I + v_cross + expm::`%^%`(v_cross, 2) * ((1-c)/(s^2))

    #  R = I + v_cross + (v_cross %*% v_cross) * ((1-c)/(s^2))
    #  R = I + v_cross + (v_cross %*% v_cross) * ((1)/(1+c))

    return(R)
  }

  ## ----------------------------------------------------------------------------
  ## Step 4: Reapply the offset
  ## Apply the transformation on the data
  R <- rotate2vertical(direction = c(a,b,c))
  #cat("\n\nRotation matrix:\n")
  #print(R)
  df_aligned <- as.data.frame(as.matrix(reduced[,1:3]) %*% R)
  ## Add the offset again
  colnames(df_aligned) <- c("x","y","z")
  df_aligned$x <- df_aligned$x + mean_x
  df_aligned$y <- df_aligned$y + mean_y
  df_aligned$z <- df_aligned$z + mean_z

  ## Store the rotation matrix for use in CloudCompare
  cc_matrix <- matrix(c(R[,1],0,R[2,],0,R[3,],0,0,0,0,1), nrow = 4, ncol = 4, byrow = T)
  write.table(cc_matrix, rotation_matrix, sep=" ", col.names = F, row.names = F)

  df_aligned_full_dataset <- as.data.frame(as.matrix(full_dataset[,1:3]) %*% R)
  colnames(df_aligned_full_dataset) <- c("x","y","z")
  df_aligned_full_dataset$x <- df_aligned_full_dataset$x #+ mean_x
  df_aligned_full_dataset$y <- df_aligned_full_dataset$y #+ mean_y
  df_aligned_full_dataset$z <- df_aligned_full_dataset$z #+ mean_z

  ## Store the transformed full dataset
  write.table(df_aligned_full_dataset, output, sep=" ", col.names = F, row.names = F)

  ## Store the transformed terrain dataset
  write.table(df_aligned, output2, sep=" ", col.names = F, row.names = F)

  ## Compute vertical angle (slope) for surface normal
  #l <- sqrt(a^2 + b^2 + c^2) # vector length, for normalization
  #an=a/l; bn=b/l; cn=c/l
  slope=acos(c)*180/pi # convert to degree
  cat("Slope:", slope,"\n")

  ## ----------------------------------------------------------------------------
  ## Step 5: Print the result

  if (plotting==T){
    s=5# scale of vector length, just for visual effect
    ## Plot original data
    plot3d(df[,1:3], col="blue", asp="iso", size=3,
           xlim=c(min(df$x)-5,max(df$x)+5),
           ylim=c(min(df$y)-5,max(df$y)+5),
           zlim=c(min(df$z)-5,max(df$z)+5))
    planes3d(a,b,c,d,alpha=0.3,col="blue")
    arrow3d(c(mean_x,mean_y,mean_z),c(a*s+mean_x,b*s+mean_y,c*s+mean_z), col="blue", type = "rotation")
    points3d(data.frame(x=mean_x,y=mean_y,z=mean_z), size=8, col="green",add=T)
    # Plot parellel data
    points3d(df_aligned, col="darkgreen", size=3, add=T)
    arrow3d(c(mean_x,mean_y,mean_z),c(mean_x+0,mean_y+0,mean_z+1*s), col="darkgreen", type = "rotation")
    planes3d(0,0,1,-mean_z,alpha=0.3,col="darkgreen")
    axes3d(box=F)
  }

  ## Remove temporary data
  rm(full_dataset)
  rm(terrain)
  gc() # clear cache/garbage

  setwd(current_wd) # Set wd back to original wd

}

