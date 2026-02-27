#' Correct perfusion and ventilation parameters based on activity levels
#'
#' @param parms Named list or vector of baseline physiological parameters (e.g., from physB)
#' @param activityReg Character vector of activity levels for each time point (e.g., "rest", "light", "moderate", "heavy")
#' @param sex Character or numeric. Sex of the individual ("male", "female", 1, or 2)
#' @param approx_fun Logical. If TRUE, returns time-interpolating functions for each parameter (default: TRUE)
#'
#' @return List of corrected parameters (either as vectors or as interpolation functions)
#' @export
f.activity <- function(parms, activityReg, sex, approx_fun=TRUE) {
  #'  Correct perfusion and ventilation parameters based on activity levels 
  #'  defined in activityReg
  #'  Input: activityReg (vector of activity levels corresponding to time points of simulation)
  #'         parms (named physB parameter vector or list)
  #'         sex (1/"male" or 2/"female")
  #'  Returns: correctedParms (corrected physB parameter list)

  # Define which parameters are perfusion/ventilation and their correction factors by name
  perfusion_names <- c(
    "Qhrt", "Qskn", "Qadp", "Qmus", "Qbon", "Qbrn", "Qthy", "Qgon",
    "Qkid", "Qsto", "Qint", "Qspl", "Qpan", "Qliv", "Qrem"
  )
  # Correction factors for each activity level (must match order above)
  factors <- list(
    rest     = rep(1, length(perfusion_names)),
    light    = c(2,1,1,5,1,1,1,1,0.75,0.75,0.75,0.75,0.75,0.75,0.75),
    moderate = c(3,1,1,10,1,1,1,1,0.5,0.5,0.5,0.5,0.5,0.5,0.5),
    heavy    = c(4,1,1,15,1,1,1,1,0.25,0.25,0.25,0.25,0.25,0.25,0.25)
  )

  correctedParms <- list()
  n_times <- length(activityReg)
  times <- seq_len(n_times)

  # Initialize Qc as zero vector
  Qc_vec <- rep(0, n_times)

  for (name in names(parms)) {
    value <- parms[[name]]
    if (name %in% perfusion_names) {
      # Get index for this perfusion parameter
      idx <- match(name, perfusion_names)
      # Build correction vector for this parameter
      corr <- numeric(n_times)
      corr[activityReg == "rest"]     <- value * factors$rest[idx]
      corr[activityReg == "light"]    <- value * factors$light[idx]
      corr[activityReg == "moderate"] <- value * factors$moderate[idx]
      corr[activityReg == "heavy"]    <- value * factors$heavy[idx]
      correctedParms[[name]] <- if (approx_fun) approxfun(times, corr, rule=2) else corr
      Qc_vec <- Qc_vec + corr
    } else {
      # Not activity-corrected: just copy
      correctedParms[[name]] <- value
    }
  }

  # Cardiac output (Qc): sum of all perfusion rates, divided by 0.95
  Qc_vec <- Qc_vec / 0.95
  correctedParms[["Qc"]] <- if (approx_fun) approxfun(times, Qc_vec, rule=2) else Qc_vec

  # Pulmonary perfusion (Qp): sex-specific formula
  if (approx_fun) {
    if (sex %in% c(1, "male")) {
      correctedParms[["Qp"]] <- approxfun(times, 0.96 * Qc_vec + 129, rule=2)
    } else {
      correctedParms[["Qp"]] <- approxfun(times, 0.62 * Qc_vec + 179, rule=2)
    }
  } else {
    if (sex %in% c(1, "male")) {
      correctedParms[["Qp"]] <- 0.96 * Qc_vec + 129
    } else {
      correctedParms[["Qp"]] <- 0.62 * Qc_vec + 179
    }
  }

  return(correctedParms)
}