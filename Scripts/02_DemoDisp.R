##### Set up functions for demography and ant dispersal ---------------------------------------------------

# Define demography function for collecting parameters into single vector
adsp.param <- function(aNum, aVal){
  
  # Prepare vector to store all demographic parameters
  aParam <- c()
  
  # Set parameters for rosette size after establishment in autumn
  aParam[1] <- demo_rose            # Size at establishment mean
  aParam[2] <- demo_rose_err        # Size at establishment SD
  
  # Set parameters for rosette survival, and growth in spring
  aParam[3] <- demo_surv            # Prob. of survival
  aParam[4] <- demo_grow_NW[1]      # Size after growth intercept
  aParam[5] <- demo_grow_NW[2]      # Size after growth size-slope
  aParam[6] <- demo_grow_err        # Size after growth SD
  
  # Set parameters for flowering probability, and flower head count
  aParam[7] <- demo_flow            # Prob. of flowering
  aParam[8] <- demo_head_NW[1]      # Head count intercept
  aParam[9] <- demo_head_NW[2]      # Head count size-slope
  aParam[10] <- demo_head_err       # Head count SD
  
  # Set parameters for flower head height above ground (in cm)
  aParam[11] <- demo_hdht_NW[1]     # Head height intercept
  aParam[12] <- demo_hdht_NW[2]     # Head height size-slope
  aParam[13] <- demo_hdht_err       # Head height SD
  
  # Set parameters for seed production, survival, establishment, and seed bank dynamics
  aParam[14] <- 476                 # Seed count per flower head
  aParam[15] <- 0.850               # Prob. of surviving pre-dispersal seed predation (florivory)
  aParam[16] <- 0.233               # Prob. of establishment from seed
  aParam[17] <- 0.233               # Prob. of seed entering seed bank if not establishing
  aParam[18] <- 0.233               # Prob. of seed establishing from seed bank
  aParam[19] <- 0.260               # Prob. of seed survival in seed bank
  aParam[20] <- 0.948               # Prob. of seed removal by ants
  aParam[21] <- 0.100               # Prob. of surviving predation if removed by ants
  
  # Scale only specified demographic parameter
  aVec <- rep(1, length(aParam))
  aVec[aNum] <- aVal
  aParam <- aParam*aVec
  
  # Return vector of demographic parameters
  return(aParam)}

# Define demography function for survival, growth, reproduction, and more
adsp.demo <- function(dType, aVec, n = 0, rsize = 0){
  
  # Import vector of demographic parameters
  aParam <- aVec
  
  # Initial rosette size upon establishment
  # Round negatives up to zero
  if(dType == "size"){
    vals <- aParam[1] + rnorm(n, mean = 0, sd = aParam[2])
    vals[vals < 0] <- 0
    return(vals)}
  
  # Rosette survival over winter
  if(dType == "surv"){
    vals <- sample(c(1, 0), size = n, prob = c(aParam[3], 1 - aParam[3]), replace = TRUE)
    return(vals)}
  
  # Rosette size after spring growth
  # Round negatives up to zero
  if(dType == "grow"){
    vals <- aParam[4] + aParam[5]*rsize + rnorm(length(rsize), mean = 0, sd = aParam[6])
    vals[vals < 0] <- 0
    return(vals)}
  
  # Rosette bolting and subsequent flowering
  if(dType == "flow"){
    vals <- sample(c(1, 0), size = n, prob = c(aParam[7], 1 - aParam[7]), replace = TRUE)
    return(vals)}
  
  # Flower head production
  # Round negatives up to zero, and any decimals to the nearest integer
  if(dType == "head"){
    vals <- exp(aParam[8] + aParam[9]*rsize + rnorm(length(rsize), mean = 0, sd = aParam[10]))
    vals[vals < 0] <- 0
    vals <- round(vals)
    return(vals)}
  
  # Flower head heights above ground
  # Round negatives up to zero, and convert height from cm to m
  if(dType == "hdht"){
    vals <- aParam[11] + aParam[12]*rsize + rnorm(length(rsize), mean = 0, sd = aParam[13])
    vals[vals < 0] <- 0
    vals <- vals/100
    return(vals)}
  
  # Per-head production of seeds, and subsequent pre-/post-dispersal seed survival
  # First term is number of viable seeds that survive florivory
  # Second term is proportion of seeds not removed by ants, or removed by ants but not eaten
  # Third term is proportion of seeds that will either establish, or not establish but enter seed bank
  # Round to nearest whole seed
  if(dType == "seed"){
    seed <- aParam[14]*aParam[15]
    prop1 <- ((1 - aParam[20]) + (aParam[20]*aParam[21]))
    prop2 <- (aParam[16] + (1 - aParam[16])*aParam[17])
    vals <- round(seed*prop1*prop2)
    return(vals)}
  
  # Establishment of seeds not entering the seed bank
  # Probability recalculated to condition on individuals that did not experience post-predation death
  # Since we already accounted for post-predation death, non-establishing seeds must enter seed bank
  if(dType == "estb"){
    prob <- (aParam[16])/(aParam[16] + (1 - aParam[16])*aParam[17])
    vals <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)
    return(vals)}
  
  # Establishment of seeds from the seed bank
  if(dType == "SBestb"){
    vals <- sample(c(1, 0), size = n, prob = c(aParam[18], 1 - aParam[18]), replace = TRUE)
    return(vals)}
  
  # Survival of seeds in the seed bank
  if(dType == "SBsurv"){
    vals <- sample(c(1, 0), size = n, prob = c(aParam[19], 1 - aParam[19]), replace = TRUE)
    return(vals)}
  
  # Removal of seeds by ants following primary dispersal
  # Probability recalculated to condition on individuals that did not experience predation
  if(dType == "ants"){
    prob <- (aParam[20]*aParam[21])/((1 - aParam[20]) + (aParam[20]*aParam[21]))
    vals <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)
    return(vals)}}

# Define function to determine whether a seed is taken to the nearest ant nest
adsp.disp <- function(d, range){
  
  # Calculate distance between seed and nests within the max range
  dist <- nestsR[abs(nestsR - d) <= range] - d
  
  # Select the closest nest (i.e. smallest absolute distance)
  dist <- dist[which.min(abs(dist))]
  if(identical(dist, numeric(0)) == TRUE){
    dist <- 0}
  
  # Return new seed location after being taken to nearest nest
  return(d + dist)}





##### Set up functions for wind dispersal -----------------------------------------------------------------

# Dispersal function for collecting wind dispersal parameters into single vector
wdsp.param <- function(wNum, wVal){
  
  # Prepare vector to store scalable dispersal parameters
  wParam <- c()
  
  # Set parameters for wind speed, seed terminal velocity, and vegetation height
  wParam[1] <- 0.056                # Prob. of seed release from capitulum
  wParam[2] <- 0.150                # Surrounding vegetation height (m)
  wParam[3] <- disp_tv[1]           # Log mean terminal velocity, lognormal dist.
  wParam[4] <- disp_tv[2]           # Log SD terminal velocity, lognormal dist.
  wParam[5] <- disp_ws[1]           # Shape wind speed, Weibull dist.
  wParam[6] <- disp_ws[2]           # Scale wind speed, Weibull dist.
  
  # Scale only specified dispersal parameter
  # Note: transformation on TV parameters transforms mean/SD, NOT log mean/SD
  # Note: if transforming wind speeds, use wNum = 5 for mean and wNum = 6 for SD
  if(wNum == 1){
    wParam[1] <- wParam[1]*wVal}
  if(wNum == 2){
    wParam[2] <- wParam[2]*wVal}
  if(wNum == 3){
    wParam[c(3, 4)] <- transform.ln(wParam[3], wParam[4], wVal, "mean")}
  if(wNum == 4){
    wParam[c(3, 4)] <- transform.ln(wParam[3], wParam[4], wVal, "sd")}
  if(wNum == 5){
    wParam[c(5, 6)] <- transform.wb(wParam[5], wParam[6], wVal, "mean")}
  if(wNum == 6){
    wParam[c(5, 6)] <- transform.wb(wParam[5], wParam[6], wVal, "sd")}
  
  # Return vector of dispersal parameters
  return(wParam)}

# Function sampling from a dispersal kernel using WALD model (Katul et al. 2005)
# H is flower head height (m), and n is number of seeds
# Code adapted from Skarpaas and Shea (2007)
wdsp.wald <- function(n, H, wVec){
  
  # Import vector of dispersal parameters
  wParam <- wVec
  
  # Get counts for seeds that are and aren't released
  nr <- sum(sample(c(1, 0), size = n, prob = c(wParam[1], 1 - wParam[1]), replace = TRUE))
  nf <- n - nr
  
  # Initialise physical constants
  K <- 0.4          # Von Karman constant
  C0 <- 3.125       # Kolmogorov constant
  
  # Initialise other fixed quantities
  Aw <- 1.3         # Ratio of sigmaw to ustar
  h <- wParam[2]    # Grass cover height
  d <- 0.7*h        # Zero-plane displacement
  z0 <- 0.1*h       # Roughness length
  zm <- 1           # Wind speed measurement height
  
  # Simulate wind dispersal if 1 or more seeds released above canopy
  if(H > h & nr > 0){
    
    # Simulate terminal velocities from lognormal distribution
    f <- rlnorm(nr, wParam[3], wParam[4])
    
    # Simulate wind speeds from Weibull distribution
    Um <- rweibull(nr, wParam[5], wParam[6])
    
    # Calculate ustar, the friction velocity
    ustar <- K*Um*(log((zm - d)/z0))^(-1)
    
    # Set up integrand for wind speed between vegetation surface and seed release height H
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
    
    # Generate inverse Gaussian distribution for seeds that are released
    # Then marginalise onto single spatial axis, assuming no dominant wind direction
    dists <- as.numeric(rinvGauss(nr, nu = nu, lambda = lambda))*cos(runif(nr, 0, 2*pi))
    
    # Assume capitulum falls in place for seeds that do not release
    # Thus, these seeds travel a distance of zero
    dists <- c(dists, rep(0, nf))
    return(dists)}
  
  # Do not simulate dispersal if no seeds are released above canopy
  if(H <= h | nr == 0){
    return(rep(0, n))}}

# Define function to estimate dispersal distances from given point
wdsp.disp <- function(n, H, wVec, d0 = 0){
  d <- wdsp.wald(n, H, wVec) + d0
  return(d)}

