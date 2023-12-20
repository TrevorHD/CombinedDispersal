##### Set up functions for distribution transformations ---------------------------------------------------

# Function to transform raw variance and/or mean, then output corresponding shape and scale
transform.wb <- function(shape, scale, fv, which.trans){
  
  # Internal functions to calculate mean and SD of Weibull distribution
  meanW <- function(shape, scale){
    return(scale*gamma(1 + 1/shape))}
  sdW <- function(shape, scale){
    return(sqrt((scale^2)*(gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)))}
  
  # Must brute force a solution since there is no easy function for inverse gamma
  
  # Create meshpoints for combinations of shape and scale
  scaleAxis <- seq(0.1, 5, by = 0.005)
  shapeAxis <- seq(0.1, 5, by = 0.005)
  scaleMesh = rep(scaleAxis, each = length(shapeAxis))
  shapeMesh = rep(shapeAxis, times = length(scaleAxis))
  
  # Let fv be the value that the mean and/or standard deviation is multiplied by
  
  # Conditional statements for the variable(s) to be transformed
  if(which.trans == "mean"){
    fm <- fv
    fs <- 1}
  if(which.trans == "sd"){
    fm <- 1
    fs <- fv}
  if(which.trans == "both"){
    fm <- fv
    fs <- fv}
  
  # Transform variables
  newMean = meanW(shape = shape, scale = scale)*fm
  newSD = sdW(shape = shape, scale = scale)*fs
  
  # Evaluate mean and SD for at meshpoints
  # Then calculate difference between estimated and supplied mean and SD
  diffs1 <- abs((mapply(meanW, shape = shapeMesh, scale = scaleMesh) - newMean)/newMean)
  diffs2 <- abs((mapply(sdW, shape = shapeMesh, scale = scaleMesh) - newSD)/newSD)
  
  # Find argmin by equally weighting mean and SD differences from supplied
  argmin <- which.min((diffs1 + diffs2)/2)
  
  # Return shape and scale at argmin
  return(c(shapeMesh[argmin], scaleMesh[argmin]))}

# Function to transform raw variance and/or mean, then output corresponding meanlog and sdlog
transform.ln <- function(meanlog, sdlog, fv, which.trans){
  
  # Internal functions to convert between mean/meanlog and sd/sdlog
  # Need to do this when transforming raw variance and/or mean
  mean.to.meanlog <- function(mean, sd){
    return(2*log(mean) - 0.5*log((mean^2) + (sd^2)))}
  sd.to.sdlog <- function(mean, sd){
    return(sqrt(-2*log(mean) + log((mean^2) + (sd^2))))}
  meanlog.to.mean <- function(meanlog, sdlog){
    return(exp(meanlog + 0.5*(sdlog^2)))}
  sdlog.to.sd <- function(meanlog, sdlog){
    return(sqrt(exp(2*meanlog + (sdlog^2))*(exp(sdlog^2) - 1)))}
  
  # Let fv be the value that the mean and/or standard deviation is multiplied by
  
  # Conditional statements for the variable(s) to be transformed
  if(which.trans == "mean"){
    fm <- fv
    fs <- 1}
  if(which.trans == "sd"){
    fm <- 1
    fs <- fv}
  if(which.trans == "both"){
    fm <- fv
    fs <- fv}
  
  # Transform variables
  newMean <- meanlog.to.mean(meanlog, sdlog)*fm
  newSD <- sdlog.to.sd(meanlog, sdlog)*fs
  newMeanlog <- mean.to.meanlog(newMean, newSD)
  newSDlog <- sd.to.sdlog(newMean, newSD)
  
  # Output new meanlog and sdlog
  return(c(newMeanlog, newSDlog))}

# Function to transform log odds ratio to probability
invlogit <- function(x){
  return(exp(x)/(1 + exp(x)))}





##### Set up functions for dispersal ----------------------------------------------------------------------

# Dispersal function for collecting parameters into single vector
wald.param <- function(sNum, sVal){
  
  # Prepare vector to store scalable dispersal parameters
  sParam <- c()
  
  # Set parameters for wind speed, seed terminal velocity, and vegetation height
  sParam[1] <- 0.15                 # Vegetation height in m
  sParam[2] <- ws_params[1]         # Shape, wind speed Weibull dist.
  sParam[3] <- ws_params[2]         # Scale, wind speed Weibull dist.
  sParam[4] <- tv_params[1]         # Log mean terminal velocity, lognormal dist.
  sParam[5] <- tv_params[2]         # Log SD terminal velocity, lognormal dist.
  sParam[6] <- an_params[1]         # Intercept, ant dispersal prob. as function of distance
  sParam[7] <- an_params[2]         # Slope, ant dispersal prob. as function of distance
  
  # Note: if transforming wind speeds, use sNum=2 for mean and sNum=3 for SD
  # Note: transformation on TV parameters transforms mean/SD, NOT log mean/SD
  if(sNum == 1){
    sParam[1] <- sParam[1]*sVal}
  if(sNum == 2){
    sParam[c(2, 3)] <- transform.wb(sParam[2], sParam[3], sVal, "mean")}
  if(sNum == 3){
    sParam[c(2, 3)] <- transform.wb(sParam[2], sParam[3], sVal, "sd")}
  if(sNum == 6){
    sParam[c(4, 5)] <- transform.ln(sParam[sNum], sParam[sNum + 1], sVal, "mean")}
  if(sNum == 7){
    sParam[c(4, 5)] <- transform.ln(sParam[sNum - 1], sParam[sNum], sVal, "sd")}
  if(sNum == 8){
    sParam[6] <- sParam[6]*sVal}
  if(sNum == 9){
    sParam[7] <- sParam[7]*sVal}
  
  # Return vector of dispersal parameters
  return(sParam)}

# Function generating a dispersal kernel using WALD model (Katul et al. 2005)
# Code adapted from Skarpaas and Shea (2007)
wald <- function(n, H, sVec){
  
  # Import vector of dispersal parameters
  sParam <- sVec
  
  # Initialise physical constants
  K <- 0.4          # von Karman constant
  C0 <- 3.125       # Kolmogorov constant
  
  # Initialise other fixed quantities
  Aw <- 1.3         # Ratio of sigmaw to ustar
  h <- sParam[1]    # Grass cover height
  d <- 0.7*h        # Zero-plane displacement
  z0 <- 0.1*h       # Roughness length
  zm <- 1           # Wind speed measurement height
  
  # Let n be the number of simulation replications
  # Let H be the seed release height
  
  # Simulate wind dispersal if released above canopy
  if(H > h){
    
    # Simulate wind speeds from Weibull distribution
    # Then marginalise onto single spatial axis, assuming no dominant wind direction
    Um <- rweibull(n, sParam[2], sParam[3])
    Um <- Um*cos(runif(n, 0, 2*pi))
    
    # Simulate terminal velocities from lognormal distribution
    f <- rlnorm(n, sParam[4], sParam[5])
    
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
    
    # Generate inverse Gaussian distribution
    # Generate more than n to deal with NAs and then cut down to n
    nInc <- ifelse(n == 1, 3, 2)
    dists <- as.numeric(na.omit(rinvGauss(n*nInc, nu = nu, lambda = lambda)))
    return(dists[1:n])}
  
  # No dispersal if released below canopy
  if(H <= h){
    return(rep(0, n))}}

# Function for probability of seed removal via ant
ant <- function(dist, sVec){
  
  # Import vector of dispersal parameters
  sParam <- sVec
  
  # Calculate probability of seed removal
  prob <- exp(sParam[6] + sParam[7]*dist)
  
  # Return probability of seed removal
  return(prob)}

# Function to see if a seed is taken to the nearest nest
nestsearch <- function(d, range, sVec){
  dists <- abs(d - nestsR)
  centre <- nestsR[which.min(dists)]
  toProb <- ant(min(dists), sVec)
  toNest <- sample(c(0, 1), 1, prob = c(1 - toProb, toProb))
  ifelse(toNest == 1 && min(dists) <= range, return(centre), return(d))}

# Function to estimate dispersal distances from given point
kern <- function(n, h, sVec, d0 = 0){
  d <- wald(n, h, sVec) + d0
  return(d)}





##### Set up functions for demography ---------------------------------------------------------------------

# Demography function for collecting parameters into single vector
demo.param <- function(dNum, dVal){
  
  # Prepare vector to store all demographic parameters
  dParam <- c()
  
  # Parameters for seed production, survival, establishment, and seed bank dynamics
  dParam[1] <- 160                  # Seeds per flower head
  dParam[2] <- 0.233                # Prob. of establishment from seed (0.302 W)
  dParam[3] <- 0.90                 # Prob. of seed predation
  dParam[4] <- 0.233                # Prob. of seed entering seed bank
  dParam[5] <- 0.233                # Prob. of seed establishing from seed bank (0.302 W)
  dParam[6] <- 0.260                # Prob. of seed survival in seed bank
  
  # Parameters for initial rosette diameter (cm) after establishment
  dParam[7] <- fits_rs[1]           # Mean rosette size
  dParam[8] <- fits_rs[2]           # SD rosette size
  
  # Parameters for flowering probability and head production as function of log rosette area
  dParam[9] <- flow_rs_NW           # Prob. flowering intercept
  dParam[10] <- 0                   # Prob. flowering slope, not applicable
  dParam[11] <- mod_head_NW[1]      # Num. heads intercept
  dParam[12] <- mod_head_NW[2]      # Num. heads slope
  
  # Parameters for distribution of flower head heights among flowering individuals
  dParam[13] <- fits_hd_NW[1]       # Mean head height
  dParam[14] <- fits_hd_NW[2]       # SD head height
  
  # Parameters for survival probability of rosettes as function of log rosette area
  dParam[15] <- surv_rs_NW          # Prob. survival intercept
  dParam[16] <- 0                   # Prob. survival slope, not applicable
  
  # Scale only specified demographic parameter
  dVec <- rep(1, length(dParam))
  dVec[dNum] <- dVal
  dParam <- dParam*dVec
  
  # Return vector of demographic parameters
  return(dParam)}

# Demography function for survival, reproduction, and more
demo <- function(dType, dVec, n = 0, rsize = 0, nflow = 0){
  
  # Import vector of demographic parameters
  dParam <- dVec
  
  # Per-head production of seeds, and subsequent seed survival
  if(dType == "seeds"){
    nseed <- dParam[1]*(1 - dParam[3])
    return(round(nseed))}
  
  # Establishment of seeds that do not enter the seed bank
  if(dType == "estAG"){
    prob <- dParam[2]
    outcomes <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)
    return(outcomes)}
  
  # Entry of seeds into the seed bank
  if(dType == "entSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(dParam[4], 1 - dParam[4]), replace = TRUE)
    return(outcomes)}
  
  # Establishment of seeds from the seed bank
  if(dType == "estSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(dParam[5], 1 - dParam[5]), replace = TRUE)
    return(outcomes)}
  
  # Survival of seeds in the seed bank
  if(dType == "surSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(dParam[6], 1 - dParam[6]), replace = TRUE)
    return(outcomes)}
  
  # Initial rosette size upon establishment
  if(dType == "rsize"){
    ros <- rtruncnorm(n, a = min(na.omit(data_rs$DM_t)),
                      mean = dParam[7], sd = dParam[8])
    return(ros)}
  
  # Flowering probability as function of initial rosette size
  if(dType == "flowering"){
    prob1 <- dParam[9] + dParam[10]*log(pi*(rsize/2)^2)
    prob0 <- 1 - prob1
    problist <- lapply(seq_len(length(prob1)), function(i) rbind(prob1, prob0)[, i])
    outcomes <- sapply(problist, sample, x = c(1, 0), size = 1, replace = TRUE)
    return(outcomes)}
  
  # Flower production as function of initial rosette size
  # Round any non-integers up to the nearest head, and negatives up to 1
  if(dType == "flowers"){
    head <- dParam[11] + dParam[12]*log(pi*(rsize/2)^2)
    head <- ifelse(head <= 0, 1, head)
    return(ceiling(head))}
  
  # Distribution of flower heights for a given individual
  # Cap min and max heights based on observational data from flower height experiment
  # Round heights to nearest cm
  if(dType == "height"){
    height <- rnorm(n, mean = dParam[13], sd = dParam[14])/100
    height[height < ht_min/100] <- ht_min/100
    height[height > ht_max/100] <- ht_max/100
    return(round(height, 2))}
  
  # Rosette survival as function of initial rosette size
  if(dType == "survival"){
    prob1 <- dParam[15] + dParam[16]*log(pi*(rsize/2)^2)
    prob0 <- 1 - prob1
    problist <- lapply(seq_len(length(prob1)), function(i) rbind(prob1, prob0)[, i])
    outcomes <- sapply(problist, sample, x = c(1, 0), size = 1, replace = TRUE)
    return(outcomes)}}

