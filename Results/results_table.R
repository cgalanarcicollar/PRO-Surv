library(xtable)

# List of filenames to load and process
files <- c("S1_p1a001.RData", "S1_p1a005.RData",  "S1_p1a01.RData",
           "S1_p05a001.RData", "S1_p05a005.RData", "S1_p05a01.RData",
           "S1_p005a001.RData", "S1_p005a005.RData", "S1_p005a01.RData",
           "S2_p1a_005.RData", "S2_p1a_01.RData",  "S2_p1a_015.RData",
           "S2_p05a_005.RData", "S2_p05a_01.RData", "S2_p05a_015.RData",
           "S2_p005a_005.RData", "S2_p005a_01.RData", "S2_p005a_015.RData")

# Create empty lists to store alpha, phi values, and final tables
alpha_phi_values_S1 <- list()
tabla_list_S1 <- list()
alpha_phi_values_S2 <- list()
tabla_list_S2 <- list()

# Loop to load files and extract the alpha, phi values
for (file in files) {
  load(file)  # Load the RData file
  dataset_name <- sub(".RData", "", file)  # Extract dataset name from the file name
  assign(dataset_name, out, envir = .GlobalEnv)  # Assign 'out' to the variable with the dataset name
  print(paste("Loaded:", dataset_name))  # Print the name of the dataset being loaded
  
  # Retrieve alpha and phi directly from the loaded dataset
  alpha <- get("alpha")  # Extract alpha
  phi <- get("phi")  # Extract phi
  
  # Store the alpha and phi values in the list with dataset name as key
  if (grepl("^S1_", dataset_name)) {
    # Store alpha and phi for S1 datasets
    alpha_phi_values_S1[[dataset_name]] <- list(alpha = alpha, phi = phi)
    
    # Add the current dataset to the S1 table list
    tabla_list_S1[[dataset_name]] <- cbind(
      get(dataset_name)$Survival.effect.table,
      100 * get(dataset_name)$prop.in.ic.alpha,
      alpha_real = alpha,  # Add alpha to the table
      phi_real = phi       # Add phi to the table
    )
  } else if (grepl("^S2_", dataset_name)) {
    # Store alpha and phi for S2 datasets
    alpha_phi_values_S2[[dataset_name]] <- list(alpha = alpha, phi = phi)
    
    # Add the current dataset to the S2 table list
    tabla_list_S2[[dataset_name]] <- cbind(
      get(dataset_name)$Survival.effect.table,
      100 * get(dataset_name)$prop.in.ic.alpha,
      alpha_real = alpha,  # Add alpha to the table
      phi_real = phi       # Add phi to the table
    )
  }
}

# Combine all the S1 datasets into one final table
tabla_S1 <- do.call(rbind, tabla_list_S1)
colnames(tabla_S1)[6] <- "CP"
# Combine all the S2 datasets into one final table
tabla_S2 <- do.call(rbind, tabla_list_S2)
colnames(tabla_S2)[6] <- "CP"



 


