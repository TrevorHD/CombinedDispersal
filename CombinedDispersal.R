##### Load libraries and initialise data ------------------------------------------------------------------

# Load libraries
library(SuppDists)
library(MASS)
library(tidyverse)
library(truncnorm)

# Load in raw weather data
data_ws1 <- read.csv("Weather1.csv")
data_ws2 <- read.csv("Weather2.csv")

# Load in terminal velocity data
# Use ambient TVs since we're only examining warming effects on height distribution
data_tv <- read.csv("SeedDropData.csv")
data_tv <- subset(data_tv, Warming == "A")





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
  return(rinvGauss(n, nu = nu, lambda = lambda))}





##### 2D expansion ----------------------------------------------------------------------------------------

plants <- data.frame(r = 0, theta = 0, x = 0, y = 0, stage = 2)

seeds <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(seeds) <- c("r", "theta", "x", "y", "germ")

# Placeholder dispersal kernel
kern <- function(n, x = 0, y = 0){
  theta <- sample(seq(0.1, 2*pi, by = pi/100), n, replace = TRUE)
  r <- rlnorm(n, meanlog = 0.25, sdlog = 1.1)
  x <- r*cos(theta) + x
  y <- r*sin(theta) + y
  return(data.frame(theta, r, x, y))}

# Placeholder nest generation
# Max range = 100 m
nests <- function(n){
  theta <- sample(seq(0.1, 2*pi, by = pi/100), n, replace = TRUE)
  r <- sample(seq(0.1, 100, by = 0.1), n, replace = TRUE)
  x <- r*cos(theta)
  y <- r*sin(theta)
  return(data.frame(theta, r, x, y))}

# for below, separate rosettes and adults first, then rejoin after performing procedures on adults

# Simulate dispersal and seed survival
for(i in 1:nrow(plants)){
  if(plants$stage[i] == 2){
    n <- demo("reproduction")
    newSeeds <- kern(n, plants$x[i], plants$y[i])
    newSeeds <- cbind(newSeeds, replicate(n, demo("germination")))
    names(newSeeds)[5] <- "germ"
    newSeeds <- newSeeds[newSeeds$germ == 1, ]
    seeds <- rbind(seeds, newSeeds)}}

# for below, separate rosettes and adults first, then rejoin after performing procedures on adults

# Kill rosettes if density is too high
for(i in 1:nrow(plants)){
  if(plants$stage[i] == 1){
    dists <- (plants$x[i] - plants$x)^2 + (plants$y[i] - plants$y)^2
    if(length(dists[dists < 1]) > 8){
      plants <- plants[-i, ]
    }
  }
}

# Adults die after reproducing
plants <- plants[plants$stage != 2, ]

# Rosettes from previous year become adults
plants$stage[plants$stage == 1] <- 2

# Surviving seeds become rosettes
seeds <- seeds[, !names(seeds) == c("germ")]
seeds$stage <- rep(1, nrow(seeds))
plants <- rbind(plants, seeds)

# Reset seeds
# Can change this later to account for seed bank
seeds <- data.frame(matrix(ncol = 5, nrow = 0))

plot(plants$x, plants$y, xlim = c(-50, 50), ylim = c(-50, 50), col = rgb(r = 0, g = 0, b = 0, alpha = 0.2), pch = 16)





##### 1D expansion ----------------------------------------------------------------------------------------

# Placeholder dispersal kernel until WALD is implemented
kern <- function(n, d0 = 0){
  d <- rlnorm(n, meanlog = 0.25, sdlog = 1.1) + d0
  return(d)}

# Initialise data; start with a single rosette and adult
plants <- data.frame(d = c(0.01, 0.01), stage = c(1, 2))
seeds <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(seeds) <- c("d", "germ")
vals <- c()

# Simulate dispersal and seed survival
for(i in 1:nrow(plants)){
  if(plants$stage[i] == 2){
    n <- demo("reproduction")
    newSeeds <- data.frame(kern(n, plants$d[i]))
    newSeeds <- cbind(newSeeds, rep(1, n))
    names(newSeeds) <- c("d", "germ")
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
  mutate(bin = as.integer(cut(d, breaks = seq(0, ceiling(max(plants$d)), 1)))) %>% 
  arrange(bin, desc(stage)) %>% 
  group_by(bin) %>%
  slice(1:pmin(10, n())) %>% 
  data.frame() -> plants
plants <- plants[, !names(plants) == c("bin")]

# Reset seeds
# Can change this later to account for seed bank
seeds <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(seeds) <- c("d", "germ")

# Plot density over space
hist(plants$d, breaks = seq(0, 500, by = 1), ylim = c(0, 12), xlab = "Distance", ylab = "Density")
vals <- c(vals, max(plants$d))

# Get wavespeeds
diff(vals)











# Demography
# In a given iteration, run reproduction, then germination, then survival, then growth
demo <- function(dType){
  if(dType == "reproduction"){
    nSeed <- 1000
    pPred <- 0.95
    return(nSeed - nSeed*pPred)}
  if(dType == "survival"){
    prob <- 0.95
    outcome <- sample(c(0, 1), 1, prob = c(1 - prob, prob))
    return(outcome)}
  if(dType == "growth"){
    return(rtruncnorm(1, a = 0.4, b = 1.7, mean = 1))}}





