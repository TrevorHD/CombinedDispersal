##### Create functions to automate elasticity calculations ------------------------------------------------

# Calculate baseline wavespeeds
wave.base <- function(){
  
  # Start system timer
  time.start <- Sys.time()
  
  # Initialise parameter vectors
  dVec1 <- demo.param(dNum = 1, dVal = 1)
  sVec1 <- wald.param(sNum = 1, sVal = 1)
  
  # Run simulation without any parameter changes
  wave1 <- wave.sim(dVec1, sVec1)
  
  # Calculate mean wavespeed
  mWave1 <- mean(diff(wave1$wavefront))
  
  # Calculate procedure time
  time.elapsed <- as.numeric(difftime(Sys.time(), time.start, units = "hours"))
  
  # Return list of calculated quantities
  return(list(procedure.time = paste0(time.elapsed, " hours"), wSpeedMean1 = mWave1,
              wSpeed1 = diff(wave1$wavefront), wFront1 = wave1$wavefront))}

# Calculate wavespeed elasticity
wave.elas <- function(dNum, dVal, sNum, sVal, baseObj){
  
  # Let dNum be the demographic parameter number
  # Let dVal be the proportion to multiply the original parameter by
  # Let baseObj be the baseline wave object
  
  # Start system timer
  time.start <- Sys.time()
  
  # Initialise parameter vectors
  if(dVal == 1 & sVal != 1){
    dVec2 <- demo.param(dNum = 1, dVal = 1)
    sVec2 <- wald.param(sNum = sNum, sVal = sVal)}
  if(dVal != 1 & sVal == 1){
    dVec2 <- demo.param(dNum = dNum, dVal = dVal)
    sVec2 <- wald.param(sNum = 1, sVal = 1)}
  
  # Run simulation with increase/decrease on a specified parameter
  wave2 <- wave.sim(dVec2, sVec2)
  
  # Calculate mean wavespeeds
  mWave1 <- baseObj$wSpeedMean1
  mWave2 <- mean(diff(wave2$wavefront))
  
  # Calculate percent change in wavespeed and parameter
  wPct <- (mWave2 - mWave1)/mWave1*100
  if(dVal == 1 & sVal != 1){
    pPct <- (sVal - 1)*100}
  if(dVal != 1 & sVal == 1){
    pPct <- (dVal - 1)*100}
  
  # Calculate elasticity
  elas <- wPct/pPct
  
  # Calculate procedure time
  time.elapsed <- as.numeric(difftime(Sys.time(), time.start, units = "hours"))
  
  # Return list of calculated quantities
  return(list(procedure.time = paste0(time.elapsed, " hours"), elasticity = elas,
              wSpeedMean1 = mWave1, wSpeedMean2 = mWave2,
              wSpeed1 = baseObj$wSpeed1, wSpeed2 = diff(wave2$wavefront),
              wFront1 = baseObj$wFront1, wFront2 = wave2$wavefront))}





##### Create function to plot wave movement ---------------------------------------------------------------

# Generate GIF of population spread
# Limit to 10000 m (nYear = 100 recommended)
generatePlots <- function(type = "hist"){
  for(i in 1:length(wave1[[2]])){
    lower <- c()
    positions <- unlist(wave1[[2]][i])
    if(length(positions) == 0){
      positions <- 0}
    if(max(positions > 500)){
      minBin <- sum(positions > floor(min(positions)) & positions < ceiling(min(positions)))
      lower <- c(rep(0:(floor(min(positions) - 1)), 10) + 0.01,
                 rep(floor(min(positions)) + 0.01, 10 - minBin))
      if(min(positions) < 1){
        lower <- sort(lower)[-c(1:20)]}}
    if(type == "hist"){
      hist(c(lower, positions), breaks = seq(0, 10000, by = 20), xlim = c(0, 10000), ylim = c(0, 250),
           xaxt = "n", yaxt = "n", xlab = "Distance (m)", ylab = "Count", main = "")
      axis(1, at = seq(0, 10000, 2000))
      axis(2, at = seq(0, 250, 50))
      text(x = 9400, y = 230, paste0("t = ", i))
      box()}
    if(type == "density"){
      plotdata <- hist(c(lower, positions), breaks = seq(0, 10000, by = 20), plot = FALSE)
      plot(plotdata$mids, plotdata$density/max(plotdata$density),
           type = "l", xlim = c(0, 10000), ylim = c(0, 1.25),
           xaxt = "n", yaxt = "n", xlab = "Distance (m)", ylab = "Relative Density", main = "")
      axis(1, at = seq(0, 10000, 2000))
      axis(2, at = seq(0, 1.25, 0.25))
      text(x = 9400, y = 1.15, paste0("t = ", i))}}}
save_gif(generatePlots("hist"), "Spread1.gif", delay = 0.3, width = 1280, height = 720, res = 144)
save_gif(generatePlots("density"), "Spread2.gif", delay = 0.3, width = 1280, height = 720, res = 144)





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
  opt_min <- ht_min_CN/100
  opt_max <- ht_max_CN/100}
if(opt_species == "CA"){
  opt_min <- ht_min_CA/100
  opt_max <- ht_max_CA/100}

# Create vector of release heights in 1 cm intervals
# Sequence bounded by minimum and maximum observed heights
opt_ht <- seq(opt_min, opt_max, by = 0.01)

# Simulate dispersal distances for each release height
for(i in 1:length(opt_ht)){
  
  # For a given release height, simulate 1 million dispersal events rounded to the nearest cm
  # Repeat 20 times, rank distances in each replicate, and average across replicates for each rank
  # This helps remove "lumpiness" of distribution, especially at higher distances
  temp <- replicate(10, round(wald(1000000, opt_ht[i], opt_species, sVec), 2))
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
    data_disp <- cbind(data_disp, temp$x)}
  
  # Remove unneeded data to free up RAM  
  remove(temp)}

# Set column names as release heights
names(data_disp) <- opt_ht

# Write full set of dispersal data to CSV
write.csv(data_disp, paste0("Data/DispersalSimulation", opt_species, ".csv"))

# Free up any unused memory
gc()

# Read dispersal data from CSV using SQL since file is large
# Standard read.csv will be too slow and may not work on machines with limited memory
data_disp <- read.csv.sql("Data/DispersalSimulationCN.csv", sql = "SELECT * FROM file")
data_disp <- data_disp[, 2:ncol(data_disp)]
names(data_disp) <- opt_ht

# Close SQL connection, and free up any unused memory
sqldf(); gc()





##### Simulate dispersal from pre-generated data ----------------------------------------------------------

# Set vector of indices to pull from dispersal heights
disp_index <- rep(1, length(opt_ht))

# Estimate dispersal distances from given point
kern <- function(n, h, sVec, disp = data_disp, d0 = 0){

# Pull column index based on height
ci <- h*100 - 15

# Pull row index for given column
index_h <- disp_index[ci]

# Pull index:(index + n) from column and save distances, then advance index by n
# Avoid use of additional functions to reduce computational load
# Also accounts for wrap-around
if(index_h + n > 1000000){
  dists <- c(disp[index_h:1000000, ci], disp[0:(index_h + n - 1000001), ci]) + d0
  disp_index[ci] <<- index_h + n - 1000000
} else {
  dists <- disp[index_h:(index_h + n - 1), ci] + d0
  disp_index[ci] <<- index_h + n}

# Return distances
return(dists)}





##### 2D expansion ----------------------------------------------------------------------------------------

# Estimate dispersal distances from given point; assume 1m plant height
kern <- function(n, x0 = 0, y0 = 0){
  r <- WALD.b(n, 1)
  theta <- sample(seq(0.1, 2*pi, by = pi/100), n, replace = TRUE)
  x <- r*cos(theta) + x0
  y <- r*sin(theta) + y0
  return(data.frame(theta, r, x, y))}

# Initialise data; start with a single rosette and adult
plants <- data.frame(r = c(0.01, 0.01), theta = c(0, 0), x = c(0.01, 0.01), y = c(0.01, 0.01), stage = c(0, 2))
seeds <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(seeds) <- c("r", "theta", "x", "y", "germ")
vals <- c()

# Placeholder nest generation
# Max range = 100 m
#nests <- function(n){
#  theta <- sample(seq(0.1, 2*pi, by = pi/100), n, replace = TRUE)
#  r <- sample(seq(0.1, 100, by = 0.1), n, replace = TRUE)
#  x <- r*cos(theta)
#  y <- r*sin(theta)
#  return(data.frame(theta, r, x, y))}

# Simulate dispersal and seed survival
for(i in 1:nrow(plants)){
  if(plants$stage[i] == 2){
    n <- demo("reproduction")
    newSeeds <- kern(n, plants$x[i], plants$y[i])
    newSeeds <- cbind(newSeeds, rep(1, n))
    names(newSeeds)[5] <- "germ"
    newSeeds <- newSeeds[newSeeds$germ == 1, ]
    seeds <- rbind(seeds, newSeeds)}}

# Kill adults after they reproduce
plants <- plants[plants$stage != 2, ]

# Rosettes from previous year become adults
plants$stage[plants$stage == 0] <- 2

# Surviving seeds become rosettes
seeds$stage <- rep(0, nrow(seeds))
seeds <- seeds[, !names(seeds) == c("germ")]
plants <- rbind(plants, seeds)

# Kill rosettes if density is too high
# Do this by sorting so that adults come first and are prioritised
plants %>% 
  mutate(xbin = as.integer(cut(x, breaks = seq(floor(min(plants$x)), ceiling(max(plants$x)), 1))),
         ybin = as.integer(cut(y, breaks = seq(floor(min(plants$y)), ceiling(max(plants$y)), 1)))) %>% 
  arrange(xbin, ybin, desc(stage)) %>% 
  group_by(xbin, ybin) %>%
  slice(1:pmin(10, n())) %>% 
  data.frame() -> plants
plants <- plants[, !names(plants) %in% c("xbin", "ybin")]

# Reset seeds
# Can change this later to account for seed bank
seeds <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(seeds) <- c("r", "theta", "x", "y", "germ")

# Plot density over space
# Note: r is relative to last point of dispersal, so we calculate distance manually
plot(plants$x, plants$y, xlim = c(-200, 200), ylim = c(-200, 200),
     col = rgb(r = 0, g = 0, b = 0, alpha = 0.2), pch = 16)
vals <- c(vals, max(sqrt((plants$x)^2 + (plants$y)^2)))

# Get wavespeeds
diff(vals)

