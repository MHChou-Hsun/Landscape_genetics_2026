library(raster)

# Import coordinates of each population
samples <- read.csv("./examples/GPS.csv")
coords <- samples[, c("Longitude", "Latitude")]

# Import environmental data
asc_files <- list.files(path = "./examples/clim_data", pattern = "*.asc", full.names = TRUE)

env_stack <- stack(asc_files)
names(env_stack)

env_values <- raster::extract(env_stack, coords)
env_data <- cbind(samples, env_values)
head(env_data)

# Assume env_data has: Populations, Latitude, Longitude, then environmental variables (from column 4 onward)
pop_names <- env_data$Populations

# Loop over each environmental variable
for (var in colnames(env_data)[4:ncol(env_data)]) {
  
  # Extract values for this variable
  values <- env_data[[var]]
  names(values) <- pop_names
  
  # Calculate pairwise Canberra distances
  dist_mat <- as.matrix(dist(values, method = "canberra"))
  
  # Keep only lower triangle
  lower_tri <- dist_mat
  lower_tri[upper.tri(lower_tri)] <- ""
  
  # Set row and column names to population names
  rownames(lower_tri) <- pop_names
  colnames(lower_tri) <- pop_names
  
  # Write CSV file
  file.name <- paste0("pairwise_", var, ".csv")
  file.DIR <- file.path("./examples/clim_data/", file.name)
  write.table(lower_tri,
              file.DIR,
              sep = ",",
              row.names = TRUE,
              col.names = NA)
}

