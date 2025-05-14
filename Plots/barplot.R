# Load required libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

# Load data
load("S1_p1a001.RData") # 'out' contains simulation results

# Convert measurement times to long format
dat.long <- pivot_longer(dat, cols = 2:5, names_to = "meas", values_to = "time")

# Combine with simulated PRO measurements
dat.y <- cbind(dat.long, y = out$y.all[, 1])
colnames(dat.y) <- c("z", "meas", "time", "y")

# Extract simulated observation times
obs_time <- out$obstime.all[, 1]

# Keep only actually observed data (based on each individual's obs time)
obs_data <- bind_rows(lapply(1:nRand, function(i) {
  dat.y[dat.y$z == i & dat.y$time <= obs_time[i], ]
}))

# Rename columns for clarity
colnames(obs_data) <- c("z", "meas", "x", "y")

# Get baseline measurement per subject (first time point)
obs_data.id <- obs_data[!duplicated(obs_data$z), ]

# Frequency table of baseline measurements
freq_table <- obs_data.id %>%
  count(y, name = "frequency") %>%
  rename(score = y)

# Plot histogram of observed measurements
histogram_plot <- obs_data %>%
  ggplot(aes(x = y)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "#e9ecef", alpha = 0.9) +
  scale_x_continuous(breaks = seq(0, 24, by = 1)) +
  theme_bw() +
  labs(x = "Measurement", y = "Frequency") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

# Display plot
histogram_plot
