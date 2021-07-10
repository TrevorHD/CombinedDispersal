##### Load libraries and initialise data ------------------------------------------------------------------

# Load libraries
library(SuppDists)
library(MASS)
library(tidyverse)
library(truncnorm)
library(xlsx)
library(gifski)
library(lme4)
library(lmerTest)

# Load in raw weather data
data_ws1 <- read.csv("Data/Weather1.csv")
data_ws2 <- read.csv("Data/Weather2.csv")

# Load in terminal velocity data
# Use ambient TVs since we're only examining warming effects on height distribution
data_tv <- read.delim("Data/SeedDropData2.txt")
data_tv <- subset(data_tv, !is.na(drop.time))

# Load height data as xlsx
data_ht <- read.xlsx("Data/ThistleData.xlsx", sheetName = "Flowers")

# Load rosette as xlsx
data_rs <- read.xlsx("Data/ThistleData.xlsx", sheetName = "General")

# Add survival indicator for rosettes that survived from transplant to summer
data_rs$Survival <- NA
data_rs$Survival[!is.na(data_rs$DM_t) & !is.na(data_rs$F)] <- 1
data_rs$Survival[!is.na(data_rs$DM_t) & is.na(data_rs$F)] <- 0

# List survival for 372 as NA since the plant was accidentally killed while cutting
data_rs$Survival[372] <- NA

# Get distribution of CA rosette sizes; is normally distributed
fits_rs_CA <- fitdistr(na.omit(subset(data_rs, Species == "CA")$DM_t), "normal")$estimate
ks.test(na.omit(subset(data_rs, Species == "CA")$DM_t),
        pnorm, mean = fits_rs_CA[1], sd = fits_rs_CA[2])

# Get distribution of CN rosette sizes; is normally distributed
fits_rs_CN <- fitdistr(na.omit(subset(data_rs, Species == "CN")$DM_t), "normal")$estimate
ks.test(na.omit(subset(data_rs, Species == "CN")$DM_t),
        pnorm, mean = fits_rs_CN[1], sd = fits_rs_CN[2])





##### Set up dispersal framework --------------------------------------------------------------------------

# Create PDF of wind speeds
# Assume no seed release occurs for wind speeds of zero, so remove zero values
ws_values <- c(data_ws1$Wind1, data_ws2$Wind1)
ws_values <- ws_values[ws_values > 0]
ws_pdf <- density(ws_values, from = min(ws_values), to = max(ws_values), bw = 0.05)
ws_mean <- mean(ws_values)

# Create PDF of CN wind speeds
# Terminal velocity is drop tube length (1.25 m) divided by drop time
tv_values_CN <- na.omit(1.25/subset(data_tv, species == "n")$drop.time)
tv_pdf_CN <- density(tv_values_CN, from = min(tv_values_CN), to = max(tv_values_CN), bw = 0.05)
tv_mean_CN <- mean(tv_values_CN)

# Create PDF of CN wind speeds
# Terminal velocity is drop tube length (1.25 m) divided by drop time
tv_values_CA <- na.omit(1.25/subset(data_tv, species == "a")$drop.time)
tv_pdf_CA <- density(tv_values_CA, from = min(tv_values_CA), to = max(tv_values_CA), bw = 0.05)
tv_mean_CA <- mean(tv_values_CA)

# Function generating a dispersal kernel using WALD model (Katul et al. 2005)
# Code adapted from Skarpaas and Shea (2007)
WALD.b <- function(n, H, species){
  
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
  if(species == "CN"){
    f <- rnorm(n, sample(tv_values_CN, size = n, replace = TRUE), tv_pdf_CN$bw)}
  if(species == "CA"){
    f <- rnorm(n, sample(tv_values_CA, size = n, replace = TRUE), tv_pdf_CA$bw)}
  
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





##### Get equations for survival and flowering ------------------------------------------------------------

# Model survival rates as function of rosette size
# Too few deaths to fit a reliable model
glmer(Survival ~ DM_t + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Survival)))
glmer(Survival ~ DM_t + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Survival)))

# Thus, rates will be estimated independent of rosette size (i.e. as a constant)
surv_rs_CA <- nrow(subset(data_rs, Species == "CA" & TRT == "NW" & Survival == 1))/
  nrow(subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Survival)))
surv_rs_CN <- nrow(subset(data_rs, Species == "CN" & TRT == "NW" & Survival == 1))/
  nrow(subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Survival)))

# Model survival rates as function of rosette size
# Too few instances of not flwoering to fit a reliable model
glmer(F ~ DM_t + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(F)))
glmer(F ~ DM_t + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(F)))

# Thus, rates will be estimated independent of rosette size (i.e. as a constant)
flow_rs_CA <- nrow(subset(data_rs, Species == "CA" & TRT == "NW" & F == 1))/
  nrow(subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(F)))
flow_rs_CN <- nrow(subset(data_rs, Species == "CN" & TRT == "NW" & F == 1))/
  nrow(subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(F)))





##### Set up demography framework -------------------------------------------------------------------------

# Demography function for growth, survival, and reproduction
demo <- function(dType, species, n = 0, rsize = 0, nflow = 0){
  
  # Fix static parameters
  seedsCA <- 83
  seedsCN <- 160
  estCA <- 0.136
  estCN <- 0.147
  sPred <- 0.95
  
  # Production of seeds, seed survival, and seed establishment
  if(dType == "seeds"){
    if(species == "CA"){
      seeds <- seedsCA*nflow*(1 - sPred)*estCA}
    if(species == "CN"){
      seeds <- seedsCN*nflow*(1 - sPred)*estCN}
    return(round(seeds))}
  
  # Initial rosette size from seed
  if(dType == "rsize"){
    if(species == "CA"){
      ros <- rtruncnorm(n, a = min(na.omit(data_rs$DM_t)),
                        mean = fits_rs_CA[1], sd = fits_rs_CA[2])}
    if(species == "CN"){
      ros <- rtruncnorm(n, a = min(na.omit(data_rs$DM_t)),
                        mean = fits_rs_CN[1], sd = fits_rs_CN[2])}
    return(ros)}
  
  # Flowering probability as function of rosette size
  if(dType == "flowering"){
    if(species == "CA"){
      prob <- flow_rs_CA}
    if(species == "CN"){
      prob <- flow_rs_CN}
    outcomes <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)
    return(outcomes)}
  
  # Flower production as function of rosette size
  if(dType == "flowers"){
    return(6 + round(rsize/100*50))}
  
  # Rosette survival as function of rosette size
  if(dType == "survival"){
    if(species == "CA"){
      prob <- surv_rs_CA}
    if(species == "CN"){
      prob <- surv_rs_CN}
    outcomes <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)
    return(outcomes)}
  
  # Height as function of rosette size
  if(dType == "growth"){
    return(0.7 + rsize/75*2)}}





##### 1D expansion ----------------------------------------------------------------------------------------

# Function to see if a seed is taken to the nearest nest
nestsearch <- function(d, range){
  dists <- abs(d - nests)
  centre <- nests[which.min(dists)]
  toNest <- sample(c(0, 1), 1, prob = c(0.05, 0.95))
  ifelse(toNest == 1 && min(dists) <= range, return(centre), return(d))}

# Estimate dispersal distances from given point; assume 1m plant height  
kern <- function(n, h, species, d0 = 0){
  d <- WALD.b(n, h, species) + d0
  return(d)}

# Choose species to model
species <- "CN"

# Set various parameters for wave model
nestOn <- TRUE    # Should ant nests be included
range <- 5        # Max detection range (m) from ant nests
nDens <- 0.1      # Ant nest density (nests/m)
trim <- TRUE      # Should core area of wave be trimmed?
trimAmt <- 500    # Distance (m) behind wavefront to trim
tDens <- 10       # Max thistle density per metre
plotOn <- FALSE   # Plot wave?

# Run invasion wave simulation
for(i in 1:100){
  
  # Initialise simulation data
  if(i == 1){
    
    # Generate nests
    nests <- sample(seq(0, 5000, by = 0.1), nDens*5000) + 0.01
    
    # Initialise data; start with a single rosette
    plants <- data.frame(d = 0.01, stage = 0, rsize = demo("rsize", species, n = 1),
                         h = 0, flow = 0, nflow = 0)
    seeds <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seeds) <- c("d", "germ")
    vals <- c()}
  
  # Surviving seeds become rosettes
  seeds$stage <- rep(0, nrow(seeds))
  seeds$rsize <- demo("rsize", species, n = nrow(seeds))
  seeds$h <- rep(0, nrow(seeds))
  seeds$flow<- rep(0, nrow(seeds))
  seeds$nflow <- rep(0, nrow(seeds))
  seeds <- seeds[, !names(seeds) == c("germ")]
  plants <- rbind(plants, seeds)
  
  # Kill rosettes if density is too high
  # Do this by sorting so that adults come first and are prioritised
  if(nrow(plants) > 0){
    plants %>% 
      mutate(bin = as.integer(cut(d, breaks = seq(0, ceiling(max(plants$d)), 1)))) %>% 
      arrange(bin, desc(stage)) %>% 
      group_by(bin) %>%
      slice(1:pmin(tDens, n())) %>% 
      data.frame() -> plants
    plants <- plants[, !names(plants) == c("bin")]}
  
  # Reset seeds
  # Can change this later to account for seed bank
  seeds <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(seeds) <- c("d", "germ")
  
  # Plot density over space
  if(plotOn == TRUE){
    if(i == 1){
      PlotList <- list()}
    PlotList[[i]] <- plants$d}
  
  # Remove rosettes that do not survive
  if(nrow(plants) > 0){
    plants <- plants[demo("survival", species, n = nrow(plants)) == 1, ]}
  
  # Trim core areas as wave progresses to save computational resources
  if(trim == TRUE & nrow(plants) > 0){
    plants <- plants[plants$d > max(plants$d) - trimAmt, ]}
  
  # Some proportion of rosettes reach adulthood in one year
  # Ones that don't will remain rosettes for another year before becoming adults
  plants$stage[plants$stage == 1] <- 2
  plants$stage[plants$stage == 0] <- sample(c(1, 2), size = length(plants$stage[plants$stage == 0]),
                                            prob = c(0.5, 0.5), replace = TRUE)
  
  # Simulate growth and flowering for adults
  plants$h[plants$stage == 2] <- demo("growth", species, rsize = plants$rsize[plants$stage == 2])
  if(sum(plants$stage == 2) > 0){
    plants$flow[plants$stage == 2] <- demo("flowering", species, n = length(plants$flow[plants$stage == 2]))}
  plants$nflow[plants$flow == 1] <- demo("flowers", species, rsize = plants$rsize[plants$flow == 1])
  
  # For flowering adults, simulate primary seed dispersal via wind
  # Also simulate seed survival and establishment
  for(j in 1:nrow(plants)){
    if(plants$flow[j] == 1){
      n <- demo("seeds", species, nflow = plants$nflow[j])
      newSeeds <- data.frame(kern(n, plants$h[j], species, plants$d[j]))
      newSeeds <- cbind(newSeeds, rep(1, n))
      names(newSeeds) <- c("d", "germ")
      newSeeds <- newSeeds[newSeeds$germ == 1, ]
      seeds <- na.omit(rbind(seeds, newSeeds))}}
  
  # Simulate secondary seed dispersal via ants
  if(nestOn == TRUE){
    seeds$d <- sapply(seeds$d, nestsearch, range = range)}
  
  # Store wavefront distance
  vals <- c(vals, max(plants$d))
  
  # Kill all adults after they reproduce
  plants <- plants[plants$stage != 2, ]}

# Get mean wavespeed
mean(diff(vals))

# Generate GIF of population spread
generatePlots <- function(type = "hist"){
  for(i in 1:length(PlotList)){
    lower <- c()
    if(length(PlotList[[i]]) == 0){
      PlotList[[i]] <- 0}
    if(max(PlotList[[i]] > 500)){
      minBin <- sum(PlotList[[i]] > floor(min(PlotList[[i]])) & PlotList[[i]] < ceiling(min(PlotList[[i]])))
      lower <- c(rep(0:(floor(min(PlotList[[i]]) - 1)), 10) + 0.01,
                 rep(floor(min(PlotList[[i]])) + 0.01, 10 - minBin))
      if(min(PlotList[[i]]) < 1){
        lower <- sort(lower)[-c(1:20)]}}
    if(type == "hist"){
      hist(c(lower, PlotList[[i]]), breaks = seq(0, 5000, by = 20), xlim = c(0, 5000), ylim = c(0, 250),
           xaxt = "n", yaxt = "n", xlab = "Distance (m)", ylab = "Count", main = "")
      axis(1, at = seq(0, 5000, 1000))
      axis(2, at = seq(0, 250, 50))
      text(x = 4700, y = 230, paste0("t = ", i))
      box()}
    if(type == "density"){
      plotdata <- hist(c(lower, PlotList[[i]]), breaks = seq(0, 5000, by = 20), plot = FALSE)
      plot(plotdata$mids, plotdata$density/max(plotdata$density), type = "l", xlim = c(0, 5000), ylim = c(0, 1.25),
           xaxt = "n", yaxt = "n", xlab = "Distance (m)", ylab = "Relative Density", main = "")
      axis(1, at = seq(0, 5000, 1000))
      axis(2, at = seq(0, 1.25, 0.25))
      text(x = 4700, y = 1.15, paste0("t = ", i))}}}
save_gif(generatePlots("hist"), "Spread1.gif", delay = 0.3, width = 1280, height = 720, res = 144)
save_gif(generatePlots("density"), "Spread2.gif", delay = 0.3, width = 1280, height = 720, res = 144)

# Extra code for 1-D nestsearch function
centres <- nests[dists < range]
dists <- dists[dists < range]
if(length(dists) > 0){
  probs <- abs(1 - ptruncnorm(dists, a = -range, b = range, mean = 0, sd = 4))*2
  for(i in 1:length(probs)){
    result <- sample(c(0, 1), 1, prob = c(probs[i], 1 - probs[i]))}}





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

