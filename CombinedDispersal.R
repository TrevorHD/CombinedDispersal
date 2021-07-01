##### Load libraries and initialise data ------------------------------------------------------------------

# Load libraries
library(SuppDists)
library(MASS)
library(tidyverse)
library(truncnorm)

# Load in raw weather data
data_ws1 <- read.csv("Data/Weather1.csv")
data_ws2 <- read.csv("Data/Weather2.csv")

# Load in terminal velocity data
# Use ambient TVs since we're only examining warming effects on height distribution
data_tv <- read.csv("Data/SeedDropData.csv")
data_tv <- subset(data_tv, Warming == "A" & Mow == "CTL")





##### Set up dispersal framework --------------------------------------------------------------------------

# Create PDF of wind speeds
# Assume no seed release occurs for wind speeds of zero, so remove zero values
ws_values <- c(data_ws1$Wind1, data_ws2$Wind1)
ws_values <- ws_values[ws_values > 0]
ws_pdf <- density(ws_values, from = min(ws_values), to = max(ws_values), bw = 0.05)
ws_mean <- mean(ws_values)

# Create PDF of wind speeds
# Terminal velocity is drop tube length divided by drop time
tv_values <- na.omit(1.25/data_tv$DT.Avg)
tv_pdf <- density(tv_values, from = min(tv_values), to = max(tv_values), bw = 0.05)
tv_mean <- mean(tv_values)

# Function generating a dispersal kernel using WALD model (Katul et al. 2005)
# Code adapted from Skarpaas and Shea (2007)
WALD.b <- function(n, H){
  
  # Initialise physical constants
  K <- 0.4      # von Karman constant
  C0 <- 3.125   # Kolmogorov constant
  
  # Initialise other fixed quantities
  Aw <- 1.3     # Ratio of sigmaw to ustar
  h <- 0.15     # Grass cover height
  d <- 0.7*h    # Zero-plane displacement
  z0 <- 0.1*h   # Roughness length
  zm <- 1       # Wind speed measurement height
  
  # Let n be the number of simulation replications
  # Let H be the seed release height
  
  # Simulate wind speeds from empirical distribution of wind speeds
  Um <- rnorm(n, sample(ws_values, size = n, replace = TRUE), ws_pdf$bw)
  
  # Simulate terminal velocities from empirical distribution of terminal velocities
  f <- rnorm(n, sample(tv_values, size = n, replace = TRUE), tv_pdf$bw)
  
  # Calculate ustar, the friction velocity
  ustar <- K*Um*(log((zm - d)/z0))^(-1)
  
  # Set up integrand for wind speed between vegetation surface and drop height H
  integrand <- function(z){
    (1/K)*(log((z - d)/z0))}
  
  # Integrate to obtain U
  U <- (ustar/H)*integrate(integrand, lower = d + z0, upper = H)$value
  
  # Calculate instability parameter
  sigma <- 2*(Aw^2)*sqrt((K*(H - d)*ustar)/(C0*U))
  
  # Calculate scale parameter lambda
  lambda <- (H/sigma)^2
  
  # Calculate location parameter nu
  nu <- H*U/f
  
  # Generate inverse Gaussian distribution
  # Generate more than n to deal with NAs and then cut down to n
  dists <- as.numeric(na.omit(rinvGauss(n*1.25, nu = nu, lambda = lambda)))
  return(dists[1:n])}





##### Set up demography framework -------------------------------------------------------------------------

# In a given iteration, run reproduction, then germination, then survival, then growth
demo <- function(dType){
  if(dType == "reproduction"){
    nSeed <- 374
    pPred <- 0.99
    return(nSeed - nSeed*pPred)}
  if(dType == "survival"){
    prob <- 0.95
    outcome <- sample(c(0, 1), 1, prob = c(1 - prob, prob))
    return(outcome)}
  if(dType == "growth"){
    return(rtruncnorm(1, a = 0.4, b = 1.7, mean = 1))}}





##### 2D expansion ----------------------------------------------------------------------------------------

# Estimate dispersal distances from given point; assume 1m plant height
kern <- function(n, x0 = 0, y0 = 0){
  r <- WALD.b(n, 1)
  theta <- sample(seq(0.1, 2*pi, by = pi/100), n, replace = TRUE)
  x <- r*cos(theta) + x0
  y <- r*sin(theta) + y0
  return(data.frame(theta, r, x, y))}

# Initialise data; start with a single rosette and adult
plants <- data.frame(r = c(0.01, 0.01), theta = c(0, 0), x = c(0.01, 0.01), y = c(0.01, 0.01), stage = c(1, 2))
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
plants$stage[plants$stage == 1] <- 2

# Surviving seeds become rosettes
seeds$stage <- rep(1, nrow(seeds))
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





##### 1D expansion ----------------------------------------------------------------------------------------

# Function to see if a seed is taken to the nearest nest
nestsearch <- function(d, range){
  dists <- abs(d - nests)
  centre <- nests[which.min(dists)]
  toNest <- sample(c(0, 1), 1, prob = c(0.05, 0.95))
  ifelse(toNest == 1 && min(dists) <= range, return(centre), return(d))}

# Estimate dispersal distances from given point; assume 1m plant height  
kern <- function(n, h, d0 = 0){
  d <- WALD.b(n, h) + d0
  return(d)}

# Generate nests; density d = 0.1 nests/m
nests <- sample(seq(0, 25000, by = 0.1), 0.1*25000)

# Initialise data; start with a single rosette
plants <- data.frame(d = 0.01, stage = 1, rsize = rnorm(1, mean = 25, sd = 3), h = 0, nflow = 0)
seeds <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(seeds) <- c("d", "germ")
vals <- c()

# Set various parameters for wave model
nestOn <- TRUE    # Should ant nests be included
range <- 5        # Max detection range (m) from ant nests
trim <- TRUE      # Should core area of wave be trimmed?
trimAmt <- 500    # Distance (m) behind wavefront to trim
tDens <- 10       # Max thistle density per metre
plotOn <- FALSE   # Plot wave?

# Run invasion wave simulation
for(i in 1:100){
  
  # Simulate primary dispersal via wind and seed survival
  for(i in 1:nrow(plants)){
    if(plants$stage[i] == 2){
      n <- round(plants$nflow[i]*demo("reproduction"))
      newSeeds <- data.frame(kern(n, plants$h[i], plants$d[i]))
      newSeeds <- cbind(newSeeds, rep(1, n))
      names(newSeeds) <- c("d", "germ")
      newSeeds <- newSeeds[newSeeds$germ == 1, ]
      seeds <- rbind(seeds, newSeeds)}}

  # Simulate secondary dispersal via ants
  if(nestOn == TRUE){
    seeds$d <- sapply(seeds$d, nestsearch, range = range)}

  # Kill adults after they reproduce
  plants <- plants[plants$stage != 2, ]

  # Rosettes from previous year become adults
  plants$h[plants$stage == 1] <- 0.25 + plants$rsize[plants$stage == 1]/100*3
  plants$nflow[plants$stage == 1] <- 6 + round(plants$rsize[plants$stage == 1]/100*50)
  plants$stage[plants$stage == 1] <- 2

  # Surviving seeds become rosettes
  seeds$stage <- rep(1, nrow(seeds))
  seeds$rsize <- rnorm(nrow(seeds), mean = 25, sd = 3)
  seeds$h <- rep(0, nrow(seeds))
  seeds$nflow <- rep(0, nrow(seeds))
  seeds <- seeds[, !names(seeds) == c("germ")]
  plants <- rbind(plants, seeds)

  # Kill rosettes if density is too high
  # Do this by sorting so that adults come first and are prioritised
  plants %>% 
    mutate(bin = as.integer(cut(d, breaks = seq(0, ceiling(max(plants$d)), 1)))) %>% 
    arrange(bin, desc(stage)) %>% 
    group_by(bin) %>%
    slice(1:pmin(tDens, n())) %>% 
    data.frame() -> plants
  plants <- plants[, !names(plants) == c("bin")]

  # Trim core areas as wave progresses to save computational resources
  if(trim == TRUE){
    plants <- plants[plants$d > max(plants$d) - trimAmt, ]}
  
  # Reset seeds
  # Can change this later to account for seed bank
  seeds <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(seeds) <- c("d", "germ")

  # Plot density over space
  if(plotOn == TRUE){
    hist(plants$d, breaks = seq(0, 1000, by = 1), ylim = c(0, 12), xlab = "Distance", ylab = "Density")}
  
  # Store wavefront distance
  vals <- c(vals, max(plants$d))}

# Get wavespeeds
diff(vals)

# Extra code for 1-D nestsearch function
centres <- nests[dists < range]
dists <- dists[dists < range]
if(length(dists) > 0){
  probs <- abs(1 - ptruncnorm(dists, a = -range, b = range, mean = 0, sd = 4))*2
  for(i in 1:length(probs)){
    result <- sample(c(0, 1), 1, prob = c(probs[i], 1 - probs[i]))}}

