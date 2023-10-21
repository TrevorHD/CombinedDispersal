##### Pre-generate dispersal data to optimise wavespeed code ----------------------------------------------

# Simulating dispersal events on-the-fly causes huge CPU bottleneck during wavespeed simulations
# However, simulating dispersal now and pulling from the data later will reduce CPU load

# Get minimum and maximum observed heights for CN and CA, from height distribution experiment
ht_min_CN <- min(subset(data_ht, Species == "CN", select = Height))
ht_max_CN <- max(subset(data_ht, Species == "CN", select = Height))
ht_min_CA <- min(subset(data_ht, Species == "CA", select = Height))
ht_max_CA <- max(subset(data_ht, Species == "CA", select = Height))

# Set baseline dispersal parameters
sVec <- wald.param(sNum = 1, sVal = 1)

# Set species before generating dispersal data; start with CN
# Need to do this manually and run all of the following code once for each species
opt_species <- "CN"
if(opt_species == "CN"){
  opt_min <- ht_min_CN
  opt_max <- ht_max_CN}
if(opt_species == "CA"){
  opt_min <- ht_min_CA
  opt_max <- ht_max_CA}

# Create vector of release heights in 1 cm intervals
# Sequence bounded by minimum and maximum observed heights
opt_ht <- seq(opt_min, opt_max, by = 0.01)

# Simulate dispersal distances for each release height
for(i in 1:length(opt_ht)){
  
  # For a given release height, simulate 1 million dispersal events
  # Repeat 20 times, rank distances in each replicate, and average across replicates for each rank
  # This helps remove "lumpiness" of distribution, especially at higher distances
  temp <- replicate(10, wald(1000000, opt_ht[i], opt_species, sVec))
  temp2 <- rowMeans(apply(temp, 2, sort, decreasing = FALSE))
  
  # Randomise order of dispersal distances
  temp2 <- sample(temp2)

  # Save distances in csv then delete R object
  # This ensures we don't run out of RAM since we're simulating so many events
  write.csv(temp2, paste0("Data/Disp_", opt_ht[i], ".csv"))
  remove(temp, temp2)}

# Combine simulated dispersal data into single dataframe
for(i in 1:length(opt_ht)){
  
  # Read dispersal distances for each release height
  temp <- read.csv(paste0("Data/Disp_", opt_ht[i], ".csv"))
  
  # Create dataframe of distances, with a column for each height
  if(i == 1){
    data_disp <- data.frame(temp$x)}
  if(i > 1){
    data_disp <- cbind(test, temp$x)}
  
  # Remove unneeded data to free up RAM  
  remove(temp)}

# Set column names as release heights
names(data_disp) <- opt_ht

# Write full set of dispersal data to CSV
write.csv(data_disp, paste0("Data/DispersalSimulation", opt_species, ".csv"))





##### Load pre-generated dispersal data -------------------------------------------------------------------

# Read dispersal data from CSV using SQL since file is large
# Standard read.csv will be too slow and may not work on machines with limited memory
data_disp <- read.csv.sql("Data/DispersalSimulation.csv", sql = "SELECT * FROM file")
names(data_disp) <- opt_ht

# Close SQL connection, and free up any unused memory
sqldf(); gc()





##### Define function for pulling from pre-generated dispersal data ---------------------------------------

# Set index for each height
disp_index <- rep(1, length(opt_ht))

# NEW: Estimate dispersal distances from given point
kern <- function(n, h, species, sVec, disp = data_disp, disp_index = disp_index, d0 = 0){
  
  # Pull column index based on height
  ci <- h*100 - 15
  
  # Pull row index for given column
  index_h <- disp_index[ci]
  
  # Pull index:(index + n) from column and save distances, then advance index by n
  # Also accounts for wrap-around
  if(index_h + n > 1000000){
    dists <- c(disp[index_h:1000000, ci], disp[0:(index_h + n - 1000001), ci]) + d0
    disp_index[ci] <- index_h + n - 1000000
  } else {
    dists <- disp[index_h:(index_h + n - 1), ci] + d0
    disp_index[ci] <- index_h + n}
  
  # Return distances and index vector as list
  return(list(dists = dists, disp_index = disp_index))}

