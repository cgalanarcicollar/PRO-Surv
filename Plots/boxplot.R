library(ggplot2)

# List of filenames to load and process
files <- c("S1_p1a001.RData", "S1_p1a005.RData",  "S1_p1a01.RData",
           "S1_p05a001.RData", "S1_p05a005.RData", "S1_p05a01.RData",
           "S1_p005a001.RData", "S1_p005a005.RData", "S1_p005a01.RData",
           "S2_p1a_005.RData", "S2_p1a_01.RData",  "S2_p1a_015.RData",
           "S2_p05a_005.RData", "S2_p05a_01.RData", "S2_p05a_015.RData",
           "S2_p005a_005.RData", "S2_p005a_01.RData", "S2_p005a_015.RData")

for (file in files) {
  load(file)  # Load the RData file
  dataset_name <- sub(".RData", "", file)  # Extract the dataset name from the file name
  assign(dataset_name, out, envir = .GlobalEnv)  # Assign 'out' to the variable with the dataset name
  print(paste("Loaded:", dataset_name))  # Print the name of the dataset being loaded
}


# Boxplot for Alpha estimation (S1 for different Phi values)
 
boxplot(S1_p005a01$TS.estimation,  S1_p005a01$TVCM.estimation, S1_p005a01$JM.estimation,
        S1_p05a01$TS.estimation,  S1_p05a01$TVCM.estimation, S1_p05a01$JM.estimation,
        S1_p1a01$TS.estimation,  S1_p1a01$TVCM.estimation, S1_p1a01$JM.estimation,
        main="Alpha estimation (S1)", ylim=c(-0.05,0.46), ylab =expression(paste(alpha,"=  0.1")),
        col=c("white","grey","grey40"))
axis(side=1, at=c(2,5,8), labels=c(expression(paste(Phi,"=0.05")),expression(paste(Phi,"=0.5")),expression(paste(Phi,"=1"))),
     line=-0.6, lwd=0, cex.axis=0.8)
abline(h=0.1, lwd=2, col="red")
legend("topleft", legend=c("TS","TVCM","JM"), bty="n", fill=c("white","grey","grey40"), border="black", pt.cex=1.75, cex=0.65)



# Boxplot for Alpha estimation (S2 for different Phi values)

boxplot(S2_p005a_01$TS.estimation,  S2_p005a_01$TVCM.estimation, S2_p005a_01$JM.estimation,
        S2_p05a_01$TS.estimation,  S2_p05a_01$TVCM.estimation, S2_p05a_01$JM.estimation,
        S2_p1a_01$TS.estimation,  S2_p1a_01$TVCM.estimation, S2_p1a_01$JM.estimation,
        main="Alpha estimation (S2)", ylim=c(-0.7,0.4), ylab =expression(paste(alpha,"= - 0.1")),
        col=c("white","grey","grey40"))
axis(side=1, at=c(2,5,8), labels=c(expression(paste(Phi,"=0.05")),expression(paste(Phi,"=0.5")),expression(paste(Phi,"=1"))),
     line=-0.6, lwd=0, cex.axis=0.8)
abline(h=-0.1, lwd=2, col="red")
legend("topleft", legend=c("TS","TVCM","JM"), bty="n", fill=c("white","grey","grey40"), border="black", pt.cex=1.75, cex=0.65)

