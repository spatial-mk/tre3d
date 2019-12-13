
#' Fit beta distribution.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Print/Display tree crown projection area.
#' @param profile A 3D point cloud file, e.g. "tree1.xyz".
#'     Assumes first three columns to be x y z without header. Futher columns are ignored.
#' @return Beta distribution
#' @examples
#' fit_beta("tree1.xyz")

## An R-Script that fits a beta distribution with optimized parameters a and b
## Matthias Kunz, 28.05.2019
## See:
## https://de.wikipedia.org/wiki/Beta-Verteilung
## http://r.789695.n4.nabble.com/How-to-fit-a-random-data-into-Beta-distribution-td3492528.html
## Beta distribution for tree profiles (from Williams et al. 2017, NEE)
## (r(a+b) / r(a)*r(b)) * x^(a-1) * (1-x)^(b-1)

## ----------------------------------------------------------------------------
## Load R packages and libraries
#library(rgl)
#library(VoxR)
#library(MASS) # fitdistr
#source("D:/Dropbox/TUD_KUNZ/tre3d/Code/tre3d/R/get_2d_chull.R")

## ----------------------------------------------------------------------------
## Sample data
# tree <- read.table("D:/Dropbox/TUD_KUNZ/R/SampleData/S10_2016_05-10.xyz",
#                  sep=" ", header = F, col.names = c("x","y","z","i"))
#
# ## Function that standardizes values between 0 and 1
# standardize_height <- function(x){(x-min(x))/(max(x)-min(x))}
#
# tree$z <- standardize_height(tree$z)
#
# ## Make voxel grid
# tree <- vox(tree, res=0.04)
# colnames(tree) <- c("x","y","z","i")
#
# ## Set lowest tree point (z axis) to zero
# tree[,3] <- tree[,3] - min(tree[,3])
#
# ## Generate strata (30 cm) based on lowest point to highest point
# slice_thickness = 0.3
# strata <- seq(0, max(tree[,3]), slice_thickness) # default 0.3 m slice
#
# ## Generate data frame that holds result for each strata
# result_slice <- data.frame(lower=strata, upper=strata + slice_thickness, radius=NA)
#
# ## Get distribution
# for (ii in result_slice$lower){
#   schicht_centroid <- tree[tree$z >= ii & tree$z < ii + slice_thickness, ]
#   schicht_centroid_xyz <- schicht_centroid[!duplicated(schicht_centroid[,1:2]),]
#   if (nrow(unique(schicht_centroid_xyz))>3 & nrow(unique(schicht_centroid_xyz[,1:2]))>3){
#       ## Convex_Hull_Variante
#       center <- get_2d_chull(unique(schicht_centroid_xyz[,1:2]), plot=F)
#       center$polygon <- NULL
#       center <- as.data.frame(center)
#   } else {
#     ## Set area and perimeter to zero if there are less than three points in slice
#     center=data.frame(area=0, perimeter=0, cpa_cen_x=pos_x, cpa_cen_y=pos_y)
#   }
#   r_ii <- sqrt(center$area/pi) ## Compute theoretical cylinder radius in slice (r_ii)
#   result_slice[result_slice$lower==ii,]$radius <- r_ii
# }
#
# result_slice$lower <- standardize_height(result_slice$lower)
# profile <- c(result_slice$radius)

## ----------------------------------------------------------------------------
## Function that return a and b from beta distribution
fit_beta <- function(profile=c(1,1.1,1)){
  plot(dist)
  lines(dist)
  beta = data.frame(a=1, b=1)
  #x<-fitdistr(dist,"beta",list(a=1,b=1))
  print(x)
  return(beta)
}

## ----------------------------------------------------------------------------
## Visualise distribution


## ----------------------------------------------------------------------------
## Test of function
#fit_beta(profile = profile)

