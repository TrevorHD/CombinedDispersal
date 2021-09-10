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





##### Set up functions for dispersal ----------------------------------------------------------------------

# Dispersal function for collecting parameters into single vector
wald.param <- function(sNum, sVal){
  
  # Prepare vector to store scalable dispersal parameters
  sParam <- c()
  
  # Set parameters for wind speed, seed terminal velocity, and vegetation height
  sParam[1] <- 0.15                 # Vegetation heignt in m
  sParam[2] <- ws_params[1]         # Shape, wind speed Weibull dist.
  sParam[3] <- ws_params[2]         # Scale, wind speed Weibull dist.
  sParam[4] <- tv_params_CA[1]      # Log mean terminal velocity, lognormal dist. (CA)
  sParam[5] <- tv_params_CA[2]      # Log SD terminal velocity, lognormal dist. (CA)
  sParam[6] <- tv_params_CN[1]      # Log mean terminal velocity, lognormal dist. (CN)
  sParam[7] <- tv_params_CN[2]      # Log SD terminal velocity, lognormal dist. (CN)
  
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
    sParam[c(6, 7)] <- transform.ln(sParam[sNum], sParam[sNum + 1], sVal, "mean")}
  if(sNum == 7){
    sParam[c(6, 7)] <- transform.ln(sParam[sNum - 1], sParam[sNum], sVal, "sd")}
  
  # Return vector of dispersal parameters
  return(sParam)}

# Function generating a dispersal kernel using WALD model (Katul et al. 2005)
# Code adapted from Skarpaas and Shea (2007)
wald <- function(n, H, species, sVec){
  
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
    
    # Simulate wind speeds from empirical distribution of wind speeds
    Um <- rweibull(n, sParam[2], sParam[3])
    
    # Simulate terminal velocities from lognormal distribution
    if(species == "CA"){
      f <- rlnorm(n, sParam[4], sParam[5])}
    if(species == "CN"){
      f <- rlnorm(n, sParam[6], sParam[7])}
    
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
    nInc <- ifelse(n == 1, 3, 2)
    dists <- as.numeric(na.omit(rinvGauss(n*nInc, nu = nu, lambda = lambda)))
    return(dists[1:n])}
  
  # No dispersal if released below canopy
  if(H <= h){
    return(rep(0, n))}}





##### Set up functions for demography ---------------------------------------------------------------------

# Demography function for collecting parameters into single vector
demo.param <- function(dNum, dVal){
  
  # Prepare vector to store all demographic parameters
  dParam <- c()
  
  # Parameters for seed production, survival, establishment, and seed bank dynamics
  dParam[1] <- 83                   # Seeds per flower head (CA)
  dParam[2] <- 160                  # Seeds per flower head (CN)
  dParam[3] <- 0.136                # Prob. of establishment from seed (CA)
  dParam[4] <- 0.147                # Prob. of establishment from seed (CN)
  dParam[5] <- 0.90                 # Prob. of seed predation
  dParam[6] <- 0.233                # Prob. of seed entering seed bank
  dParam[7] <- 0.233                # Prob. of seed establishing from seed bank
  dParam[8] <- 0.260                # Prob. of seed survival in seed bank
  
  # Parameters for initial rosette sizes (cm) after establishment
  dParam[9] <- fits_rs_CA[1]        # Mean rosette size (CA)
  dParam[10] <- fits_rs_CA[2]       # SD rosette size (CA)
  dParam[11] <- fits_rs_CN[1]       # Mean rosette size (CN)
  dParam[12] <- fits_rs_CN[2]       # SD rosette size (CN)
  
  # Parameters for flowering probability and head production as function of rosette size
  dParam[13] <- flow_rs_CA          # Prob. flowering intercept (CA)
  dParam[14] <- 0                   # Prob. flowering slope, not applicable (CA)
  dParam[15] <- flow_rs_CN          # Prob. flowering intercept (CN)
  dParam[16] <- 0                   # Prob. flowering slope, not applicable (CN)
  dParam[17] <- mod_head_CA         # Num. heads intercept (CA)
  dParam[18] <- 0                   # Num. heads slope, not applicable (CA)
  dParam[19] <- mod_head_CN[1]      # Num. heads intercept (CN)
  dParam[20] <- mod_head_CN[2]      # Num. heads slope (CN)
  
  # Parameters for distribution of flower head heights among flowering individuals
  dParam[21] <- fits_hd_CA_NW[1]    # Mean head height (CA)
  dParam[22] <- fits_hd_CA_NW[2]    # SD head height (CA)
  dParam[23] <- fits_hd_CN_NW[1]    # Mean head height (CN)
  dParam[24] <- fits_hd_CN_NW[2]    # SD head height (CA)
  
  # Parameters for survival probability of rosettes as function of rosette size
  dParam[25] <- surv_rs_CA          # Prob. survival intercept (CA)
  dParam[26] <- 0                   # Prob. survival slope, not applicable (CA)
  dParam[27] <- surv_rs_CN          # Prob. survival intercept (CN)
  dParam[28] <- 0                   # Prob. survival slope, not applicable (CN)
  
  # Scale only specified demographic parameter
  dVec <- rep(1, length(dParam))
  dVec[dNum] <- dVal
  dParam <- dParam*dVec
  
  # Return vector of demographic parameters
  return(dParam)}

# Demography function for survival, reproduction, and more
demo <- function(dType, species, dVec, n = 0, rsize = 0, nflow = 0){
  
  # Import vector of demographic parameters
  dParam <- dVec
  
  # Per-head production of seeds, and subsequent seed survival
  if(dType == "seeds"){
    if(species == "CA"){
      nseed <- dParam[1]*(1 - dParam[5])}
    if(species == "CN"){
      nseed <- dParam[2]*(1 - dParam[5])}
    return(round(nseed))}
  
  # Establishment of seeds that do not enter the seed bank
  if(dType == "estAG"){
    if(species == "CA"){
      prob <- dParam[3]}
    if(species == "CN"){
      prob <- dParam[4]}
    outcomes <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)
    return(outcomes)}
  
  # Entry of seeds into the seed bank
  if(dType == "entSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(dParam[6], 1 - dParam[6]), replace = TRUE)
    return(outcomes)}
  
  # Establishment of seeds from the seed bank
  if(dType == "estSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(dParam[7], 1 - dParam[7]), replace = TRUE)
    return(outcomes)}
  
  # Survival of seeds in the seed bank
  if(dType == "surSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(dParam[8], 1 - dParam[8]), replace = TRUE)
    return(outcomes)}
  
  # Initial rosette size upon establishment
  if(dType == "rsize"){
    if(species == "CA"){
      ros <- rtruncnorm(n, a = min(na.omit(data_rs$DM_t)),
                        mean = dParam[9], sd = dParam[10])}
    if(species == "CN"){
      ros <- rtruncnorm(n, a = min(na.omit(data_rs$DM_t)),
                        mean = dParam[11], sd = dParam[12])}
    return(ros)}
  
  # Flowering probability as function of initial rosette size
  if(dType == "flowering"){
    if(species == "CA"){
      prob1 <- dParam[13] + dParam[14]*rsize}
    if(species == "CN"){
      prob1 <- dParam[15] + dParam[16]*rsize}
    prob0 <- 1 - prob1
    problist <- lapply(seq_len(length(prob1)), function(i) rbind(prob1, prob0)[,i])
    outcomes <- sapply(problist, sample, x = c(1, 0), size = 1, replace = TRUE)
    return(outcomes)}
  
  # Flower production as function of initial rosette size
  # Round any non-integers up to the nearest head
  if(dType == "flowers"){
    if(species == "CA"){
      head <- dParam[17] + dParam[18]*rsize}
    if(species == "CN"){
      head <- dParam[19] + dParam[20]*rsize}
    return(ceiling(head))}
  
  # Distribution of flower heights for a given individual
  if(dType == "height"){
    if(species == "CA"){
      height <- rnorm(n, mean = dParam[21], sd = dParam[22])/100}
    if(species == "CN"){
      height <- rnorm(n, mean = dParam[23], sd = dParam[24])/100}
    return(height)}
  
  # Rosette survival as function of initial rosette size
  if(dType == "survival"){
    if(species == "CA"){
      prob1 <- dParam[25] + dParam[26]*rsize}
    if(species == "CN"){
      prob1 <- dParam[27] + dParam[28]*rsize}
    prob0 <- 1 - prob1
    problist <- lapply(seq_len(length(prob1)), function(i) rbind(prob1, prob0)[,i])
    outcomes <- sapply(problist, sample, x = c(1, 0), size = 1, replace = TRUE)
    return(outcomes)}}

