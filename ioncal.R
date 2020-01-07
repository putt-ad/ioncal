
#' High Pressure Ion Chromatography Calibration Function
#'
#' USAGE INSTRUCTIONS
#' Read .csv files from Chromeleon software
#' to start clean up the header and the footer of the file.
#' Follow the style guide available on the calibration_style_guide.txt



ioncal <- function(df) {
  if (!is.data.frame (df))
    warning("Did not provide a usable data frame")
  oas <- readline(prompt="name of calibration standard: ")
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
  y<- as.numeric(df_oas[[2]])

  # build a bestfit linear model for each analyte and saving each coef to `a` print summary for each
  for (i in 4:ncol(df_oas)){
    m <- lm(y~0 + as.numeric(df_oas[[i]]), na.action=na.omit)
    a[i-3] <- m$coefficients[1]
    r2[i-3] <- summary(m)$r.squared
    print(summary(m))
  }

  # Creating matrix of unknown values that need calibrated
  df_unknown <- oa_dataset[which(df$Type == "Unknown"), ]
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
  # replace calculated concentrations with measured values
  df_unknown[ ,c(4:ncol(df_unknown))] <-Y

  # Generating data output
  datalist <- list(col_names, a, r2)
  overview <- data.frame(matrix(unlist(datalist), nrow=length(datalist), byrow=TRUE))
  overview <- cbind(X0 = c("analyte", "coef", "r.squared"), overview)
  colnames(overview) <- NULL
  write.csv(df_unknown, "calibrated_analytes.csv")
  write.csv(overview, "stats_overview.csv")

  #transpose overview
  #turn datframe into matrix
  jpeg(filename = "calibration_summary.jpeg", width = 6, height = 6, units = 'in',
       pointsize = 12, quality = 75, bg = "white", res = 300)
  layout(matrix(c(1,2), 2, 2, byrow = TRUE),
         widths=c(1,1), heights=c(1,2))
  boxplot(r2, ylab = "R Squared")
  boxplot(a, ylab = "Coefficients")
  mtext("Calibration Summary", side=3, outer=TRUE, line=-3)
  dev.off()
  cat("\nSuccess! \nObtained calibration standard concentrations: ", y, "\nFor analytes:", col_names, "\nYou have produced the following output files: \ncalibrated data output: calibrated_analytes.csv \ncalbiration summary figure: calibration_summary.jpeg\nstats overview: stats_overview.csv\n\n")
}

