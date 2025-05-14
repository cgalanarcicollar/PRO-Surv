library(xtable)

# List of RDS files to process
rds_files <- list.files("Results", pattern = "^S2_.*\\.rds$", full.names = TRUE)

# Initialize list to store final tables
summary_tables <- list()

# Helper function to compute coverage probability
proportion <- function(b, se, true, level = .95, df = Inf){
  qtile <- level + (1 - level)/2
  lower.bound <- b - qt(qtile, df = df)*se
  upper.bound <- b + qt(qtile, df = df)*se
  in.ci <- ifelse(true >= lower.bound & true <= upper.bound, 1, 0)
  cp <- mean(in.ci)
  return(cp * 100)
}

# Loop through each RDS file
for (file in rds_files) {
  out <- readRDS(file)
  dataset_name <- gsub("\\.rds$", "", basename(file))
  alpha <- out$parameters[[1]]
  
  # JM estimates
  Alpha_JM <- apply(out$result_JM$surv$Alpha, 2, mean)
  asd_alphJM <- Alpha_JM[3]
  esd_alphJM <- sd(out$result_JM$surv$Alpha[, 1])
  bias_alphJM <- if (alpha == 0) Alpha_JM[1] - alpha else (Alpha_JM[1] - alpha)/alpha
  CP_JM <- mean(out$result_JM$surv$Alpha[,4] < alpha & alpha < out$result_JM$surv$Alpha[,8]) * 100
  
  # TS estimates
  Alpha_TS <- apply(out$result_TS$survival$Alpha, 2, mean, na.rm = TRUE)
  asd_alphTS <- Alpha_TS[2]
  esd_alphTS <- sd(out$result_TS$survival$Alpha[,1], na.rm = TRUE)
  bias_alphTS <- if (alpha == 0) Alpha_TS[1] - alpha else (Alpha_TS[1] - alpha)/alpha
  CP_TS <- proportion(out$result_TS$survival$Alpha[,1], out$result_TS$survival$Alpha[,2], alpha)
  
  # Create result table
  result_table <- rbind(
    JM = c(bias_alphJM, esd_alphJM, asd_alphJM, CP_JM),
    TS = c(bias_alphTS, esd_alphTS, asd_alphTS, CP_TS)
  )
  
  colnames(result_table) <- c("Bias", "ESD", "ASD", "CP")
  summary_tables[[dataset_name]] <- result_table
}

# Combine all results into one big table with identifiers
final_table <- do.call(rbind, lapply(names(summary_tables), function(name) {
  cbind(Dataset = name, Method = rownames(summary_tables[[name]]), summary_tables[[name]])
}))

# Print the final summary table
print(xtable(final_table, digits = c(0, 0, 0, 3, 3, 3, 1)), type = "html")
