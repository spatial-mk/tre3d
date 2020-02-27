## A function that prepares a FARO .fls file or .xyz for SSCI computation
## Note: Plane adjustment for scans on sloped terrain is included (optional)
## Matthias Kunz, 15.10.2019

#' SSCI computation after Ehbrecht et al. 2017
#' @author Matthias Kunz, last updated: 18.12.2019
#' @description Compute SSCI for .fls, .ptx or .xyz files
#' @param input_fls A .fls, .ptx or .xyz file. Note .fls requires a Main file in the same folder as the .fls file.
#' @param align_to_terrain Logical. Should a slope correction for SSCI be performed. Default TRUE.
#' @param voxel_size Voxel size (in meter) for ENL computation. Default 0.1 m.
#' @param slice_thickness Thickness (in meter) of vertical layers for ENL computation. Default 0.5 m.
#' @param angular_resolution Angular resolution of the scan, i.e. 360/number of points. Default 360/2560.
#' @param CROP_SOR_XYZ Logical. Should a statistical outlier removal filter be applied and point cloud be cropped. Default FALSE.
#' @param crop_square_length Half edge length (in meter) of squared crop along xy axis around scan center. Default +-10.
#' @param crop_circular Logical. Should points be clipped circular around the (0,0,0) position. Default FALSE.
#' @param crop_circular_radius Clip radius (in meter) if crop_circular=TRUE
#' @param save_steps Logical. Should aligned point clouds and filter clouds be saved. Default TRUE.
#' @param plotting Logical. Should the MeanFrac polygons be plotted. Note: This may takes a long time. Default FALSE.
#' @return List with SSCI metrics.
#' @export
#' @examples
#' SSCI()

## Required R libraries and additional packages
#library(stringr) # grep strings
#library(readr) # read files fast
#library(dplyr) # filter and select from data frame
#library(crayon) # Fancy cat output
#source("D:/Dropbox/TUD_KUNZ/R/SSCI/align_slope_plane_fit.R")
#source("D:/Dropbox/TUD_KUNZ/R/SSCI/MeanFrac.R")
#source("D:/Dropbox/TUD_KUNZ/R/SSCI/ENL.R")

# read multiple files at onces: library(data.table); df<-rbindlist(lapply(files, fread))

#if (!dir.exists("C:/LAStools/")){stop("Please install LAStools under C:/LAStools")}

## ****************************************************************************
## Function that converts and computes SSCI from .xyz/.ptx/.fls files
SSCI <- function(input_fls="*.fls", align_to_terrain=T,
                                     voxel_size=0.1, slice_thickness=0.5,
                                     angular_resolution=360/2560,CROP_SOR_xyz=F,crop_square_length=10,
                                     crop_circular=F,crop_circular_radius=10,
                                     save_steps=F,
                                     plotting=F){

  cat(crayon::red("\n--------------------"),crayon::bold(crayon::blue("\n>>> INPUT:", input_fls,"\n")))
  if (file.exists(input_fls)==F){cat(input_fls," does not exist.");return(0)}

  ## --------------------------------------------------------------------------
  ## Step 0: Do some settings and checks

  ## Check if input is a .fls/.ptx or .xyz file (PTX/FLS are treated the same)
  extension <- tools::file_ext(basename(input_fls))
  case_fls=F # Default input is .xyz
  if (extension=="fls" | extension=="ptx"){case_fls=T}
  else if (extension=="xyz"){case_fls=F}
  else {stop("\nFunction needs .fls/.ptx or .xyz input")}

  ## Set some file paths and get the name of the file
  cc <- "C://Progra~1//CloudCompare//CloudCompare.exe"  ## CloudCompare should be installed under C:/Program Files/CloudCompare
  cur_wd<-getwd() # Current working directoy
  setwd(dirname(input_fls)) ## Set working directory to .fls file location
  fls_file <- tools::file_path_sans_ext(basename(input_fls))## Get the name of the .fls/.xyz file

  ## Check if SSCI was already computed
  SSCI_out <- paste0(getwd(),"/", tools::file_path_sans_ext(basename(input_fls)),"_SSCI.txt")
  if (file.exists(SSCI_out)){cat(SSCI_out,"already computed\n"); setwd(cur_wd); return()}

  ## Helper function for circular crop
  circleFun <- function(center = c(0,0),radius = crop_circular_radius, npoints = 72){
    r = radius
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  circle_poly <- tidyr::unite(round(circleFun(),6),col="pairs",sep=" ")
  circle_poly <- paste0(circle_poly$pairs, collapse = " ")

  ## ----------------------------------------------------------------------------
  ## FLS/PTX file processing
  if (case_fls==T){

    ## Pre-check if file was already processed
    if (file.exists(paste0(getwd(),"/",fls_file, "_CROPPED_SOR.xyz"))==F){

      ## --------------------------------------------------------------------------
      ## Step 1: Convert .fls to ASCII (.bin) and crop using CloudCompare
      ## Run CloudCompare in silent mode using cmd
      cat(crayon::bold(crayon::blue("\n>>> Convert .fls/.ptx to ASCII (.xyz) and crop\n")))
      if (crop_circular==F){
        system("cmd.exe", input = paste0(cc, " -SILENT",
                                         " -O ", paste0(input_fls),
                                         " -AUTO_SAVE OFF", " -NO_TIMESTAMP",
                                         " -CROP ",-crop_square_length,":",-crop_square_length,":-1000:",crop_square_length,":",crop_square_length,":1000",
                                         #" -CROP -10:-10:-1000:10:10:1000",
                                         " -C_EXPORT_FMT",  " ASC", " -PREC 6", " -EXT", " xyz", " -SAVE_CLOUDS",
                                         " -LOG_FILE ", paste0(getwd(),"/",fls_file, "_log_crop.txt")))
      } else {
        system("cmd.exe", input = paste0(cc, " -SILENT",
                                         " -O ", paste0(input_fls),
                                         " -AUTO_SAVE OFF", " -NO_TIMESTAMP",
                                         " -CROP2D Z 72 ", circle_poly,
                                         " -C_EXPORT_FMT",  " ASC", " -PREC 6", " -EXT", " xyz", " -SAVE_CLOUDS",
                                         " -LOG_FILE ", paste0(getwd(),"/",fls_file, "_log_crop.txt")))
      }

      ## Wait until CloudCompare is finished and file CROPPED.xyz is available
      Sys.sleep(10) # Wait until CloudCompare is initialized
      check = "running"
      while(file.exists(paste0(getwd(),"/",fls_file, "_log_crop.txt"))==FALSE){
        cat("\nCloudCompare is not finished yet. Wait 60 sec longer. ", format(Sys.time(),usetz = TRUE))
        Sys.sleep(60)
      }
      unlink(paste0(getwd(),"/",fls_file, "_log_crop.txt"))

      ## --------------------------------------------------------------------------
      ## Step 2: Adjust point cloud for scanner elevation so that scan center is at 0,0,0
      ## Get the scanner elevation from the transformation information, 4x4 rotation matrix
      ## Note: Only works with newer .fls files, not with .fls from FARO Photon where center is already 0,0,0
      if (extension=="fls"){
        main_file <- list.files(path = dirname(input_fls), pattern = "Main")
        lines <- readLines(main_file, warn = F)
        id <- grep("Transformation", lines)
        transformation_info <- lines[id]
        transformation_info <- gsub("\\s+", " ", transformation_info)
        scanner_elevation <- unlist(strsplit(transformation_info, split = " "))[16] # Grep the 16th item [3,3], contains the height of the sensor
        scanner_elevation <- as.numeric(scanner_elevation)
      } else {
        scanner_elevation = 0 # No barometric height stored in FARO Photon .fls, thus height is 0
      }

      cat("\nRead in converted and cropped file\n")
      temp <- suppressMessages(readr::read_delim(paste0(fls_file, "_CROPPED.xyz"), delim = " ", col_names = F))

      ## Apply the elevation offset
      cat(crayon::bold(crayon::blue("\n>>> Adjust scanner elevation of",round(scanner_elevation,2),"m to 0 m\n")))
      temp$X3 <- temp$X3 - scanner_elevation

      ## Remove extreme points 60 m above lowest point, this sometimes happens
      temp <- dplyr::filter(temp, X3 < 50)

      ## Remove the possible RGB columns to reduce the size
      temp$X4 <- NULL; temp$X5 <- NULL; temp$X6 <- NULL

      ## Save the reduced point cloud
      write.table(x = temp, file = paste0(fls_file, "_CROPPED.xyz"), sep = " ", quote = F, col.names = F, row.names = F)
      rm(temp)

      ## --------------------------------------------------------------------------
      ## Step 3: Apply stastical outlier removal filter (SOR)
      ## #"C:", "&& cd ", cc, "&& CloudCompare.exe"
      cat(crayon::bold(crayon::blue("\n>>> Apply SOR filter\n")))
      system("cmd.exe", input = paste0(cc, " -SILENT" ,
                                       " -O ",  paste0(getwd(),"/",fls_file, "_CROPPED.xyz"),
                                       " -AUTO_SAVE ", " OFF ", " -NO_TIMESTAMP",
                                       " -SOR 10 3 ",
                                       " -C_EXPORT_FMT ", " ASC", " -PREC 6", " -EXT", " xyz ", " -SAVE_CLOUDS",
                                       " -LOG_FILE ", paste0(getwd(),"/",fls_file, "_log_sor.txt")))

      ## Wait until CloudCompare is finished and file _SOR.xyz is available
      Sys.sleep(10) # Wait until CloudCompare is initialized
      check = "running"
      while(file.exists(paste0(getwd(),"/",fls_file, "_log_sor.txt"))==FALSE){
        cat("\nCloudCompare is not finished yet. Wait 60 sec longer. ", format(Sys.time(),usetz = TRUE))
        Sys.sleep(60)
      }
      unlink(paste0(getwd(),"/",fls_file, "_log_sor.txt"))

    } # End-of-fls file processing

  } # end of FLS/PTX case if

  ## ----------------------------------------------------------------------------
  ## Step 4: Align to terrain/slope (optional)

  ## Substep: Extract terrain points using LAStools
  if (extension=="fls" | extension=="ptx"){
    cat(crayon::bold(crayon::blue("\n>>> Convert xyz to las\n")))
    system("cmd.exe", input = paste0("C:\\LAStools\\bin\\txt2las.exe ","-i ",
                                     paste0(getwd(),"/",fls_file, "_CROPPED_SOR.xyz"), " -quiet", " -parse xyz -set_scale 0.001 0.001 0.001",
                                     " -o ", paste0(getwd(),"/",fls_file, ".laz")))
  } else {

    ## Apply CROP and SOR filter xyz file if CROP_SOR_xyz == T
    if (CROP_SOR_xyz == T){
      ## CROP ---------------
      cat(crayon::bold(crayon::blue("\n>>> Crop xyz\n")))
      if (crop_circular==F){
        system("cmd.exe", input = paste0(cc, " -SILENT",
                                         " -O ", paste0(input_fls),
                                         " -AUTO_SAVE OFF", " -NO_TIMESTAMP",
                                         " -CROP ",-crop_square_length,":",-crop_square_length,":-1000:",crop_square_length,":",crop_square_length,":1000",
                                         " -C_EXPORT_FMT",  " ASC", " -PREC 6", " -EXT", " xyz", " -SAVE_CLOUDS",
                                         " -LOG_FILE ", paste0(getwd(),"/",fls_file, "_log_crop.txt")))
      } else {
        system("cmd.exe", input = paste0(cc, " -SILENT",
                                         " -O ", paste0(input_fls),
                                         " -AUTO_SAVE OFF", " -NO_TIMESTAMP",
                                         " -CROP2D Z 72 ", circle_poly,
                                         " -C_EXPORT_FMT",  " ASC", " -PREC 6", " -EXT", " xyz", " -SAVE_CLOUDS",
                                         " -LOG_FILE ", paste0(getwd(),"/",fls_file, "_log_crop.txt")))
      }

      ## Wait until CloudCompare is finished and file CROPPED.xyz is available
      Sys.sleep(10) # Wait until CloudCompare is initialized
      check = "running"
      while(file.exists(paste0(getwd(),"/",fls_file, "_log_crop.txt"))==FALSE){
        cat("\nCloudCompare is not finished yet. Wait 60 sec longer. ", format(Sys.time(),usetz = TRUE))
        Sys.sleep(60)
      }
      unlink(paste0(getwd(),"/",fls_file, "_log_crop.txt"))

      ## SOR filter ---------
      cat(crayon::bold(crayon::blue("\n>>> Apply SOR filter\n")))
      system("cmd.exe", input = paste0(cc, " -SILENT" ,
                                       " -O ",  paste0(getwd(),"/",fls_file, "_CROPPED.xyz"),
                                       " -AUTO_SAVE ", " OFF ", " -NO_TIMESTAMP",
                                       " -SOR 10 3 ",
                                       " -C_EXPORT_FMT ", " ASC", " -PREC 6", " -EXT", " xyz ", " -SAVE_CLOUDS",
                                       " -LOG_FILE ", paste0(getwd(),"/",fls_file, "_log_sor.txt")))

      ## Wait until CloudCompare is finished and file _SOR.xyz is available
      Sys.sleep(10) # Wait until CloudCompare is initialized
      check = "running"
      while(file.exists(paste0(getwd(),"/",fls_file, "_log_sor.txt"))==FALSE){
        cat("\nCloudCompare is not finished yet. Wait 60 sec longer. ", format(Sys.time(),usetz = TRUE))
        Sys.sleep(60)
      }
      unlink(paste0(getwd(),"/",fls_file, "_log_sor.txt"))

      ## CONVERT ------------
      cat(crayon::bold(crayon::blue("\n>>> Convert xyz to las\n")))
      system("cmd.exe", input = paste0("C:\\LAStools\\bin\\txt2las.exe ","-i ",
                                       paste0(getwd(),"/",fls_file, "_CROPPED_SOR.xyz"), " -quiet", " -parse xyz -set_scale 0.001 0.001 0.001",
                                       " -o ", paste0(getwd(),"/",fls_file, ".laz")))


    } else {
      cat(crayon::bold(crayon::blue("\n>>> Convert xyz to las\n")))
      system("cmd.exe", input = paste0("C:\\LAStools\\bin\\txt2las.exe ","-i ",
                                       input_fls, " -quiet", " -parse xyz -set_scale 0.001 0.001 0.001",
                                       " -o ", paste0(getwd(),"/",fls_file, ".laz")))
    }


  }
  unlink(paste0(getwd(),"/",fls_file, "_CROPPED.xyz")) ## Remove the temporary file

  ## Classify ground (and vegetation) points
  cat(crayon::bold(crayon::blue("\n>>> Classify ground and vegetation\n")))
  system("cmd.exe", input = paste0("C:\\LAStools\\bin\\lasground.exe ","-i ",
                                   paste0(getwd(),"/",fls_file, ".laz"), " -o ", paste0(getwd(),"/",fls_file, "_classification.laz"),
                                   " -quiet -all_returns -not_airborne -v -step 1 -spike 0.3"))

  ## Extract the ground points
  cat(crayon::bold(crayon::blue("\n>>> Extract terrain/ground points\n")))
  system("cmd.exe", input = paste0("C:\\LAStools\\bin\\lasground.exe ","-i ",
                                   paste0(getwd(),"/",fls_file, "_classification.laz"),
                                   " -keep_xy -5 -5 5 5", " -quiet" ," -keep_class 2 -oparse xyz",
                                   " -o ", paste0(getwd(),"/",fls_file, "_ground.laz")))
  unlink(paste0(getwd(),"/",fls_file, "_classification.laz")) ## Remove the temporary files

  ## Reduce point number before ground extraction to avoid lastools licence issues
  cat(crayon::bold(crayon::blue("\n>>> Thin point cloud to avoid LAStools licence issues\n")))
  system("cmd.exe", input = paste0("C:\\LAStools\\bin\\lasthin.exe ","-i ",
                                   paste0(getwd(),"/",fls_file, "_ground.laz"),
                                   " -step 0.04 -quiet -central -o ",
                                   paste0(getwd(),"/",fls_file, "_ground_reduced.laz")))

  ## Extract the ground points as DEM
  system("cmd.exe", input = paste0("C:\\LAStools\\bin\\las2dem.exe ","-i ",
                                   paste0(getwd(),"/",fls_file, "_ground_reduced.laz"), " -step 0.2", " -o ",
                                   paste0(getwd(),"/",fls_file, "_ground_dem.xyz")))
  terrain<-read.table(paste0(getwd(),"/",fls_file, "_ground_dem.xyz"), sep=",", header = F)## Save terrain file with space as delimeter
  write.table(x = terrain[complete.cases(terrain),], file = paste0(getwd(),"/",fls_file, "_ground_dem.xyz"), sep=" ", col.names = F, row.names = F, quote = F)

  unlink(paste0(getwd(),"/",fls_file, "_ground.laz")) ## Remove the temporary files
  unlink(paste0(getwd(),"/",fls_file, "_ground_reduced.laz")) ## Remove the temporary files
  unlink(paste0(getwd(),"/",fls_file, ".laz")) ## Remove the temporary files

  ## ----------------------------------------------------------------------------
  ## Step 5: Make alignement to terrain

  # Make alignment
  cat(crayon::bold(crayon::blue("\n>>> Align point cloud to terrain slope\n")))
  if (extension=="fls" | extension=="ptx"){
    slope<-tre3d::align_slope_plane_fit(terrain_data = paste0(getwd(),"/",fls_file, "_ground_dem.xyz"),
                        full_data = paste0(getwd(),"/",fls_file, "_CROPPED_SOR.xyz"), plotting=plotting)
  } else if (CROP_SOR_xyz==T){
    slope<-tre3d::align_slope_plane_fit(terrain_data = paste0(getwd(),"/",fls_file, "_ground_dem.xyz"),
                          full_data = paste0(getwd(),"/",fls_file, "_CROPPED_SOR.xyz"), plotting=plotting)
  } else {
    slope<-tre3d::align_slope_plane_fit(terrain_data = paste0(getwd(),"/",fls_file, "_ground_dem.xyz"),
                          full_data = input_fls, plotting=plotting)
  }

  ## ----------------------------------------------------------------------------
  ## Step 6: Compute MeanFrac, ENL and SSCI for aligned and not-aligned version

  ## Original case
  cat(crayon::bold(crayon::blue("\n>>> Compute MeanFrac, ENL and SSCI for not-aligned (original) version\n")))
  if (extension=="fls"| extension=="ptx"){
    cloud<-suppressMessages(read_delim(paste0(getwd(),"/",fls_file, "_CROPPED_SOR.xyz"),
                                       delim = " ", col_names = F))
    cloud_file <- paste0(getwd(),"/",fls_file, "_CROPPED_SOR.xyz")
  } else if (CROP_SOR_xyz==T){
    cloud<-suppressMessages(read_delim(paste0(getwd(),"/",fls_file, "_CROPPED_SOR.xyz"),
                                       delim = " ", col_names = F))
    cloud_file <- paste0(getwd(),"/",fls_file, "_CROPPED_SOR.xyz")
  } else {
    cloud<-suppressMessages(read_delim(input_fls,
                                       delim = " ", col_names = F))
    cloud_file <- input_fls
  }

  rMeanFrac <- MeanFrac(cloud, plotting = plotting, angular_resolution = angular_resolution)
  rENL <- ENL(input_cloud = cloud,  voxel_size = voxel_size, slice_thickness = slice_thickness)

  ## Compute height normalized vertical layers
  rENL_HN <- tre3d::ENL_height_normalized(input_cloud = cloud_file, voxel_size = voxel_size, slice_thickness = slice_thickness)

  SSCI <- rMeanFrac$result_MeanFrac^(log(as.numeric(rENL$ENL)))
  rm(cloud)
  gc() # Garbage (cache) collection

  ## Aligned case
  cat(crayon::bold(crayon::blue("\n>>> Compute MeanFrac, ENL and SSCI for aligned version\n")))
  if (extension=="fls"| extension=="ptx"){
    cloud_aligned<-suppressMessages(read_delim(paste0(getwd(),"/",fls_file, "_CROPPED_SOR_aligned.xyz"),
                                               delim = " ", col_names = F))
  } else if (CROP_SOR_xyz==T) {
    cloud_aligned<-suppressMessages(read_delim(paste0(getwd(),"/",fls_file, "_CROPPED_SOR_aligned.xyz"),
                                               delim = " ", col_names = F))
  } else {
    cloud_aligned<-suppressMessages(read_delim(paste0(tools::file_path_sans_ext(basename(input_fls)),"_aligned.xyz"),
                                               delim = " ", col_names = F))
  }

  rMeanFrac_aligned <- tre3d::MeanFrac(cloud_aligned, plotting = plotting, angular_resolution = angular_resolution)
  rENL_aligned <- tre3d::ENL(input_cloud = cloud_aligned, voxel_size = voxel_size, slice_thickness = slice_thickness)
  SSCI_aligned <- rMeanFrac_aligned$result_MeanFrac^(log(as.numeric(rENL_aligned$ENL)))
  rm(cloud_aligned)
  gc() # Garbage (cache) collection

  ## Delete steps (temporary file outputs) if set
  if(save_steps==F){
    unlink(paste0(getwd(),"/",fls_file, "_CROPPED_SOR.xyz"))
    unlink(paste0(tools::file_path_sans_ext(basename(input_fls)),"_aligned.xyz"))
    unlink(paste0(getwd(),"/",fls_file, "_CROPPED_SOR_aligned.xyz"))
    unlink(paste0(getwd(),"/",fls_file, "_ground_dem.xyz"))
    unlink(paste0(getwd(),"/",fls_file, "_ground_dem_aligned.xyz"))
    unlink(paste0(getwd(),"/",fls_file, "_cc_matrix.txt"))
  }

  ## Store a file with the result
  result<-data.frame(file=basename(input_fls),
                     MeanFrac=rMeanFrac$result_MeanFrac, ENL=rENL$ENL, SSCI=SSCI,
                     MeanFrac_aligned=rMeanFrac_aligned$result_MeanFrac, ENL_aligned=rENL_aligned$ENL, SSCI_aligned=SSCI_aligned,
                     voxel_size=rENL$voxel_size, slice_thickness=rENL$slice_thickness, Slope=slope, ENL_HN=rENL_HN$ENL)
  write.table(result, file = SSCI_out,
              sep = ";", col.names = T, row.names = F, quote = F)

  setwd(cur_wd) ## Set back to old working directory

  return(list(File=basename(input_fls),MeanFrac=rMeanFrac, MeanFrac_aligned=rMeanFrac_aligned,
              ENL=rENL, ENL_aligned=rENL_aligned,
              SSCI=SSCI, SSCI_aligned=SSCI_aligned, Slope=slope, ENL_HN=rENL_HN))

} ## END-OF-FUNCTION
