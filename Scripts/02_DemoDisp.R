##### Set up functions for wind dispersal -----------------------------------------------------------------

# Dispersal function for collecting wind dispersal parameters into single vector
wdsp.param <- function(wNum, wVal){
  
  # Prepare vector to store scalable dispersal parameters
  wParam <- c()
  
  # Set parameters for wind speed, seed terminal velocity, and vegetation height
  wParam[1] <- 0.15                 # Vegetation height (m)
  wParam[2] <- disp_ws[1]           # Shape wind speed, Weibull dist.
  wParam[3] <- disp_ws[2]           # Scale wind speed, Weibull dist.
  wParam[4] <- disp_tv[1]           # Log mean terminal velocity, lognormal dist.
  wParam[5] <- disp_tv[2]           # Log SD terminal velocity, lognormal dist.
  wParam[6] <- 0.056                # Probability of seed release from capitulum
  
  # Note: if transforming wind speeds, use wNum=2 for mean and wNum=3 for SD
  # Note: transformation on TV parameters transforms mean/SD, NOT log mean/SD
  if(wNum == 1){
    wParam[1] <- wParam[1]*wVal}
  if(wNum == 2){
    wParam[c(2, 3)] <- transform.wb(wParam[2], wParam[3], wVal, "mean")}
  if(wNum == 3){
    wParam[c(2, 3)] <- transform.wb(wParam[2], wParam[3], wVal, "sd")}
  if(wNum == 4){
    wParam[c(4, 5)] <- transform.ln(wParam[wNum], wParam[wNum + 1], wVal, "mean")}
  if(wNum == 5){
    wParam[c(4, 5)] <- transform.ln(wParam[wNum - 1], wParam[wNum], wVal, "sd")}
  if(wNum == 8){
    wParam[6] <- wParam[6]*wVal}
  
  # Return vector of dispersal parameters
  return(wParam)}

# Function sampling from a dispersal kernel using WALD model (Katul et al. 2005)
# H is flower head height (m), and n is number of seeds
# Code adapted from Skarpaas and Shea (2007)
wdsp.wald <- function(n, H, wVec){
  
  # Import vector of dispersal parameters
  wParam <- wVec
  
  # Initialise physical constants
  K <- 0.4          # von Karman constant
  C0 <- 3.125       # Kolmogorov constant
  
  # Initialise other fixed quantities
  Aw <- 1.3         # Ratio of sigmaw to ustar
  h <- wParam[1]    # Grass cover height
  d <- 0.7*h        # Zero-plane displacement
  z0 <- 0.1*h       # Roughness length
  zm <- 1           # Wind speed measurement height
  
  # Get counts for seeds that are and aren't released
  nr <- round(n*wParam[6], 0)
  nf <- n - nr
  
  # Simulate wind dispersal if released above canopy
  if(H > h){
    
    # Simulate wind speeds from Weibull distribution
    Um <- rweibull(nr, wParam[2], wParam[3])
    
    # Simulate terminal velocities from lognormal distribution
    f <- rlnorm(nr, wParam[4], wParam[5])
    
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

# Function to estimate dispersal distances from given point
wdsp.disp <- function(n, h, wVec, d0 = 0){
  d <- wdsp.wald(n, h, wVec) + d0
  return(d)}





##### Set up functions for demography and ant dispersal ---------------------------------------------------

# Demography function for collecting parameters into single vector
adsp.param <- function(aNum, aVal){
  
  # Prepare vector to store all demographic parameters
  aParam <- c()
  
  # Parameters for seed production, survival, establishment, and seed bank dynamics
  aParam[1] <- 476                  # Seeds per flower head
  aParam[2] <- 0.850                # Prob. of surviving pre-dispersal seed predation (florivory)
  aParam[3] <- 0.960                # Prob. of seed removal by ants
  aParam[4] <- 0.100                # Prob. of surviving post-dispersal seed predation by ants
  aParam[5] <- 0.233                # Prob. of establishment from seed
  aParam[6] <- 0.233                # Prob. of seed entering seed bank
  aParam[7] <- 0.233                # Prob. of seed establishing from seed bank
  aParam[8] <- 0.260                # Prob. of seed survival in seed bank
  
  # Parameters for log initial rosette area after establishment
  aParam[9] <- demo_rose            # Rosette size t0 intercept
  aParam[10] <- demo_rose_err       # Rosette size t0 model SD
  
  # Parameters for rosette survival, and growth as function of log rosette area
  aParam[11] <- demo_surv           # Prob. survival
  aParam[12] <- demo_grow_NW[1]     # Rosette size t1 intercept
  aParam[13] <- demo_grow_NW[2]     # Rosette size t1 size-slope
  aParam[14] <- demo_grow_err       # Mean growth model SD
  
  # Parameters for flowering probability, and head production as function of log rosette area
  aParam[15] <- demo_flow           # Prob. flowering
  aParam[16] <- demo_head_NW[1]     # Num. heads intercept
  aParam[17] <- demo_head_NW[2]     # Num. heads size-slope
  aParam[18] <- demo_head_err       # Num. heads model SD
  
  # Parameters for flower head height as function of rosette diameter
  aParam[19] <- demo_hdht_NW[1]     # Head height intercept
  aParam[20] <- demo_hdht_NW[2]     # Head height diameter-slope
  aParam[21] <- demo_hdht_err       # Head height model SD
  
  # Scale only specified demographic parameter
  aVec <- rep(1, length(aParam))
  aVec[aNum] <- aVal
  aParam <- aParam*aVec
  
  # Return vector of demographic parameters
  return(aParam)}

# Demography function for survival, reproduction, and more
adsp.demo <- function(dType, aVec, n = 0, rsize = 0){
  
  # Import vector of demographic parameters
  aParam <- aVec
  
  # Per-head production of seeds, and subsequent pre-/post-dispersal seed survival
  # Product of first two terms represents number of viable seeds that survive florivory
  # Third term represents proportion of seeds not removed by ants or removed but not eaten
  # Fourth term represents proportion of seeds that will either establish or enter seed bank
  # Round up to nearest whole seed
  if(dType == "seeds"){
    nseed <- aParam[1]*aParam[2]*((1 - aParam[3]) + (aParam[3]*aParam[4]))*(1 - aParam[5] - aParam[6])
    return(ceiling(nseed))}
  
  # Removal of seeds by ants following primary dispersal
  # Probability recalculated to condition on individuals that did not experience predation
  if(dType == "antRem"){
    prob1 <- (aParam[3]*aParam[4])/((1 - aParam[3]) + (aParam[3]*aParam[4]))
    outcomes <- sample(c(1, 0), size = n, prob = c(prob1, 1 - prob1), replace = TRUE)
    return(outcomes)}
  
  # Establishment of seeds not entering the seed bank
  # Probability recalculated to condition on individuals that did not experience post-predation death
  # We already accounted for post-predation death, so non-establishing seeds must enter seed bank
  if(dType == "estNew"){
    prob1 <- (aParam[5])/(aParam[5] + aParam[6])
    outcomes <- sample(c(1, 0), size = n, prob = c(prob1, 1 - prob1), replace = TRUE)
    return(outcomes)}
  
  # Establishment of seeds from the seed bank
  if(dType == "estSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(aParam[7], 1 - aParam[7]), replace = TRUE)
    return(outcomes)}
  
  # Survival of seeds in the seed bank
  if(dType == "surSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(aParam[8], 1 - aParam[8]), replace = TRUE)
    return(outcomes)}
  
  # Initial rosette size upon establishment
  if(dType == "rsize"){
    ros <- aParam[9] + rnorm(n, mean = 0, sd = aParam[10])
    return(ros)}
  
  # Rosette survival before growth stage
  if(dType == "surv"){
    outcomes <- sample(c(1, 0), size = n, prob = c(aParam[11], 1 - aParam[11]), replace = TRUE)
    return(outcomes)}
  
  # Rosette growth from t0 to t1
  if(dType == "grow"){
    size <- aParam[12] + aParam[13]*rsize + rnorm(length(rsize), mean = 0, sd = aParam[14])
    return(size)}
  
  # Flowering probability as function of initial rosette size
  if(dType == "flow"){
    outcomes <- sample(c(1, 0), size = n, prob = c(aParam[15], 1 - aParam[15]), replace = TRUE)
    return(outcomes)}
  
  # Flower production as function of initial rosette size
  # Round any non-integers up to the nearest head, and negatives up to 1
  if(dType == "heads"){
    head <- aParam[16] + aParam[17]*rsize + rnorm(length(rsize), mean = 0, sd = aParam[18])
    head <- ifelse(head <= 0, 1, head)
    return(ceiling(head))}
  
  # Distribution of flower heights for a given individual
  # Convert height from cm to m before returning values
  if(dType == "height"){
    height <- aParam[19] + aParam[20]*area.i(rsize) + rnorm(length(rsize), mean = 0, sd = aParam[21])
    return(height/100)}}

# Function to see if a seed is taken to the nearest ant nest
adsp.disp <- function(d, range){
  
  # Calculate distance between seed and nests within the max range
  dist <- nestsR[abs(nestsR - d) <= range] - d
  
  # Select the closest nest (i.e. smallest absolute distance)
  dist <- dist[which.min(abs(dist))]
  
  if(identical(dist, numeric(0)) == TRUE){
    dist <- 0}
  
  # Return new seed location after being taken to nearest nest
  return(d + dist)}

