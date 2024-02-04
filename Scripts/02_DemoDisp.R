##### Set up functions for dispersal ----------------------------------------------------------------------

# Dispersal function for collecting parameters into single vector
wald.param <- function(sNum, sVal){
  
  # Prepare vector to store scalable dispersal parameters
  sParam <- c()
  
  # Set parameters for wind speed, seed terminal velocity, and vegetation height
  sParam[1] <- 0.15                 # Vegetation height in m
  sParam[2] <- ws_params[1]         # Shape wind speed, Weibull dist.
  sParam[3] <- ws_params[2]         # Scale wind speed, Weibull dist.
  sParam[4] <- tv_params[1]         # Log mean terminal velocity, lognormal dist.
  sParam[5] <- tv_params[2]         # Log SD terminal velocity, lognormal dist.
  sParam[6] <- an_params[1]         # Intercept, ant dispersal prob. as function of distance
  sParam[7] <- an_params[2]         # Slope, ant dispersal prob. as function of distance
  sParam[8] <- 0.056                # Probability of seed release from capitulum
  
  # Note: if transforming wind speeds, use sNum=2 for mean and sNum=3 for SD
  # Note: transformation on TV parameters transforms mean/SD, NOT log mean/SD
  if(sNum == 1){
    sParam[1] <- sParam[1]*sVal}
  if(sNum == 2){
    sParam[c(2, 3)] <- transform.wb(sParam[2], sParam[3], sVal, "mean")}
  if(sNum == 3){
    sParam[c(2, 3)] <- transform.wb(sParam[2], sParam[3], sVal, "sd")}
  if(sNum == 4){
    sParam[c(4, 5)] <- transform.ln(sParam[sNum], sParam[sNum + 1], sVal, "mean")}
  if(sNum == 5){
    sParam[c(4, 5)] <- transform.ln(sParam[sNum - 1], sParam[sNum], sVal, "sd")}
  if(sNum == 6){
    sParam[6] <- sParam[6]*sVal}
  if(sNum == 7){
    sParam[7] <- sParam[7]*sVal
  if(sNum == 8){
    sParam[8] <- sParam[8]*sVal}}
  
  # Return vector of dispersal parameters
  return(sParam)}

# Function generating a dispersal kernel using WALD model (Katul et al. 2005)
# H is flower head height (m), and n is number of seeds
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
  
  # Get counts for seeds that are and aren't released
  nr <- round(n*sParam[8], 0)
  nf <- n - nr
  
  # Simulate wind dispersal if released above canopy
  if(H > h){
    
    # Simulate wind speeds from Weibull distribution
    Um <- rweibull(nr, sParam[2], sParam[3])
    
    # Simulate terminal velocities from lognormal distribution
    f <- rlnorm(nr, sParam[4], sParam[5])
    
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
    
    # Generate inverse Gaussian distribution for seeds that release
    # Then marginalise onto single spatial axis, assuming no dominant wind direction
    dists <- as.numeric(rinvGauss(nr, nu = nu, lambda = lambda))*cos(runif(nr, 0, 2*pi))
    
    # For seeds that do not release, assume capitilum falls in place
    # Thus, these seeds travel a distance of zero
    dists <- c(dists, rep(0, nf))
    return(dists)}
  
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
  dParam[1] <- 476                  # Seeds per flower head
  dParam[2] <- 0.15 + (0.85*0.90)   # Prob. of seed predation (pre- plus post-dispersal)
  dParam[3] <- 0.233                # Prob. of establishment from seed
  dParam[4] <- 0.233                # Prob. of seed entering seed bank
  dParam[5] <- 0.233                # Prob. of seed establishing from seed bank
  dParam[6] <- 0.260                # Prob. of seed survival in seed bank
  
  # Parameters for log initial rosette area after establishment
  dParam[7] <- mod_rose             # Rosette size t0 intercept
  dParam[8] <- mod_rose_err         # Rosette size t0 SD
  
  # Parameters for rosette survival, and growth as function of log rosette area
  dParam[9] <- surv_rs              # Prob. survival
  dParam[10] <- mod_grow_NW[1]      # Rosette size t1 intercept
  dParam[11] <- mod_grow_NW[2]      # Rosette size t1 SD
  dParam[12] <- mod_grow_err        # Mean growth SD
  
  # Parameters for flowering probability, and head production as function of log rosette area
  dParam[13] <- flow_rs             # Prob. flowering
  dParam[14] <- mod_head_NW[1]      # Num. heads intercept
  dParam[15] <- mod_head_NW[2]      # Num. heads size-slope
  dParam[16] <- mod_head_err        # Num. heads SD
  
  # Parameters for flower head height as function of rosette diameter
  dParam[17] <- mod_hdht_NW[1]      # Head height intercept
  dParam[18] <- mod_hdht_NW[2]      # Head height diameter-slope
  dParam[19] <- mod_hdht_err        # Head height SD
  
  # Scale only specified demographic parameter
  dVec <- rep(1, length(dParam))
  dVec[dNum] <- dVal
  dParam <- dParam*dVec
  
  # Return vector of demographic parameters
  return(dParam)}

# Demography function for survival, reproduction, and more
demo <- function(dType, dVec, n = 0, rsize = 0){
  
  # Import vector of demographic parameters
  dParam <- dVec
  
  # Per-head production of seeds, and subsequent seed survival
  if(dType == "seeds"){
    nseed <- dParam[1]*(1 - dParam[2])
    return(ceiling(nseed))}
  
  # Establishment of seeds that do not enter the seed bank
  if(dType == "estAG"){
    outcomes <- sample(c(1, 0), size = n, prob = c(dParam[3], 1 - dParam[3]), replace = TRUE)
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
    ros <- dParam[7] + rnorm(n, mean = 0, sd = dParam[8])
    return(ros)}
  
  # Rosette survival before growth stage
  if(dType == "survival"){
    outcomes <- sample(c(1, 0), size = n, prob = c(dParam[9], 1 - dParam[9]), replace = TRUE)
    return(outcomes)}
  
  # Rosette growth from t0 to t1
  if(dType == "grow"){
    size <- dParam[10] + dParam[11]*rsize + rnorm(length(rsize), mean = 0, sd = dParam[12])
    return(size)}
  
  # Flowering probability as function of initial rosette size
  if(dType == "flowering"){
    outcomes <- sample(c(1, 0), size = n, prob = c(dParam[13], 1 - dParam[13]), replace = TRUE)
    return(outcomes)}
  
  # Flower production as function of initial rosette size
  # Round any non-integers up to the nearest head, and negatives up to 1
  if(dType == "flowers"){
    head <- dParam[14] + dParam[15]*rsize + rnorm(length(rsize), mean = 0, sd = dParam[16])
    head <- ifelse(head <= 0, 1, head)
    return(ceiling(head))}
  
  # Distribution of flower heights for a given individual
  if(dType == "height"){
    height <- dParam[17] + dParam[18]*area.i(rsize) + rnorm(length(rsize), mean = 0, sd = dParam[19])
    return(round(height, 2)/100)}}

