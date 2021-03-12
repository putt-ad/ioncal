#' High Performance Ion Chromatography Calibration Function'
#' USAGE INSTRUCTIONS'
#' Read .csv files from Chromeleon software'
#' to start clean up the header and the footer of the in-put file.'
#' Follow the style guide available on the calibration_style_guide.txt'

ioncal <- function(df) {
  if (!is.data.frame (df))
    warning("Did not provide a usable data frame")
  file_name <- readline(prompt="name of your dataset: ")
  oas <- readline(prompt="name of calibration standard: ")
  Unknown <- readline(prompt="name of Unknown Samples: ")
  units <- readline(prompt="Concentration Units: ")
  # select operation: subset rows based on type "oas", and keep all columns
  df_oas <- df[which(df$Type == oas), ]
  #view(df_oas)
  # extracting and saving names of analytes
  col_names <- names(df_oas)[c(4:ncol(df_oas))]
  #view(col_names)
  # initialize list of coef
  a <- rep(0,length(col_names))
  # initializing list of r2 values
  r2 <- rep(0, length(col_names))
  #initializing list of standard concentrations
  y <- as.numeric(df_oas[[2]])
  
  # set pdf for calibration output
  pdf(paste0(file_name,"_calibration curves.pdf"))
  # set up a grid for plotting
  par(mfrow = c(2,2))
  # build a bestfit linear model for each analyte calibration and saving each coef to `a` print summary for each
  for (i in 4:ncol(df_oas)){
    m <- lm(y~0 + as.numeric(df_oas[[i]]), na.action=na.omit)
    a[i-3] <- m$coefficients[1]
    r2[i-3] <- summary(m)$r.squared
    print(summary(m))
    plot(df_oas[[i]] ~y, df_oas, main = colnames(df_oas)[i], xlab = paste0("Standard Conc.(",units,")"), ylab = paste0("Measured Amount(",units,")"))
    abline(lm(df_oas[[i]] ~y, df_oas))
  }
  dev.off()
  
  # Creating matrix of unknown values that need calibrated
  df_unknown <- df[which(df$Type == Unknown), ]
  # make a matrix of all analytes
  X <- data.matrix(df_unknown[c(4:ncol(df_unknown))])
  # replace NA with 0
  X[is.na(X)] <- 0
  # diagnal matrix constructed of calculated coef
  coef_diag <-diag(a, length(a), length(a))
  # do the math to turn detection into concentrations
  Y <- X %*% coef_diag
  # reassign NA
  Y[Y == 0] <- NA
  # replace measured values with calculated concentrations
  df_unknown[ ,c(4:ncol(df_unknown))] <-Y
  
  # Creating a matrix of the calibration data to be plotted
  X_cal <- data.matrix(df_oas[c(4:ncol(df_oas))])
  # replace NA with 0
  X_cal[is.na(X_cal)] <- 0
  # diag matrix for 
  coef_diag_cal <-diag(a, length(a), length(a))
  # do the math to turn detection into concentrations
  Y_cal <- X_cal %*% coef_diag_cal
  # replace the zeroes with NA
  Y_cal[Y_cal == 0] <- NA
  # replace uncalibrated data with 
  df_oas[ ,c(4:ncol(df_oas))] <-Y_cal
  
  # Generating data output
  datalist <- list(col_names, a, r2)
  overview <- data.frame(matrix(unlist(datalist), nrow=length(datalist), byrow=TRUE))
  overview <- cbind(X0 = c("analyte", "coef", "r.squared"), overview)
  colnames(overview) <- NULL
  # get working directory
  wd <- getwd()
  write.csv(df_unknown, paste0(file_name, "_calibrated_analytes.csv"))
  write.csv(overview, paste0(file_name, "_stats_overview.csv"))
  write.csv(df_oas, paste0(file_name, "_standard_results.csv"))
  #result message
  cat("\nSuccess! \nObtained calibration standard concentrations: ", y, 
      "\nFor analytes:", col_names, "\n\nThe following output files have been added to your directory(",wd,")
      \ncalibrated data output:",file_name,"calibrated_analytes.csv 
      \ncalbiration plots:",file_name,"calibration curves.pdf
      \nstats overview:" ,file_name,"stats_overview.csv\nCalculated standards:" ,file_name,"standard_results.csv")
}


