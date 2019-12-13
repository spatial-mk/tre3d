## ****************************************************************************
## Function 'CBHfromQSM' that extracts crown base height (CBH) from a QSM
## Matthias Kunz, 09.04.2019
##
## Requires: cyl_data_*.txt that contains QSM (from TreeQSM 2.30)
## ****************************************************************************

#' Crown base height (CBH) extraction from a quantitative structure model (QSM).
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Extracts the CBH from a QSM file. QSM file (cyl_data*.txt) must come from TreeQSM 2.30.
#' @references \url{https://github.com/InverseTampere/TreeQSM}
#' @param cyl_data A file containing a QSM, e.g. cyl_data_*.txt
#' @param min_diameter Minimum diameter (in meter) of a first order branch considered to be first living trunk. Default 0.07
#' @param max_angle_deg Maximum allowed branch angle (in degree) from vertical tree axis. Default 60.0
#' @param min_branch_length_m Minimum length of a first order branch considered to be first living trunk. Default 2.0
#' @return CBH in meter. Note that CBH computation sets lowest tree point in z axis to 0.
#' @export
#' @examples
#' cbh_from_qsm("cyl_data_*.txt", 0.07, 60.0, 2.0)

## ----------------------------------------------------------------------------
## CBHfromQSM function
cbh_from_qsm <- function(cyl_data="cyl_data_*.txt",
                       min_diameter_m=0.07, max_angle_deg=60,
                       min_branch_length_m=2.0){

  ## ----------------------------------
  ## Check if file exists. If not, give warning
  if (!file.exists(cyl_data)){
    stop(paste0("QSM file: ", cyl_data, " does not exist."))
  }

  ## Read the QSM file, with Try-Catch-Read-File
  mtry <- try(read.table(cyl_data, sep=",", header = T), silent = TRUE)

  if (class(mtry) != "try-error") {
    ## Lese QSM ein
    qsm <- read.table(cyl_data, sep=",", header = T)

    ## Get tree name
    tree <- tools::file_path_sans_ext(qsm)
    tree <- gsub("cyl_data_", "", tree)

    ## Set lowest point to zero to get a correct CBH
    qsm$start_3 <- qsm$start_3 - min(qsm[qsm$BranchOrder==0,]$start_3)

    ## Get angle (measured from vertical axis) of cylinders
    qsm$z_angle = rad2deg(atan2(sqrt(qsm$axis_1^2 + qsm$axis_2^2), qsm$axis_3))

    ## Select first cylinders of first branch order cylinders
    candidates <- dplyr::filter(qsm, PositionInBranch==1, BranchOrder==1)

    ## Select branches that fulfil minimum radius condition
    candidates <- dplyr::filter(candidates, radius >= min_diameter_m/2)

    ## Select branches with angle lower than 60 degree from z axis
    candidates <- dplyr::filter(candidates, z_angle < max_angle_deg)

    ## Check if there are candidates
    if (nrow(candidates)==0){
      CBH=0
    } else {

      ## Get length of candidate branches
      candidate_branches <- unique(candidates$branch)
      length_candidates <- doBy::summaryBy(length ~ branch, data = qsm[qsm$branch %in% candidate_branches,], FUN = sum)

      ## Select only those with specific length
      length_candidates <- dplyr::filter(length_candidates, length.sum > min_branch_length_m)
      candidates <- dplyr::filter(candidates, branch %in% length_candidates$branch)

      ## Select lowest start_3 (Z-Value) of candidates
      crown_base <- dplyr::select(dplyr::filter(candidates, start_3 == min(start_3)), start_1, start_2, start_3, radius)

      ## If no suitable candidate height, set CBH to zero
      if (nrow(crown_base)==0){
        CBH <- 0
      } else {
        # Get crown base height (CBH)
        CBH <- crown_base$start_3
      }

    } # End of candidate check



  } else {
    message(paste0(cyl_data," cannot be read. Please check."))
  }

  ## Return the CBH value
  return(CBH)

}

## EOF ##
