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

# List survival for 372 as NA since the plant was accidentally killed while trimming
data_rs$Survival[372] <- NA

# Get number of flowers per plant
data_ht %>% 
  subset(Type == "f" | Type == "s") %>% 
  group_by(Row, Group, Plant, TRT) %>% 
  summarise(Heads = n()) %>% 
  data.frame() -> heads
data_rs <- merge(data_rs, heads, by = c("Row", "Group", "Plant"), all = TRUE)
data_rs <- subset(data_rs, !is.na(Species), select = -c(TRT.y))
names(data_rs)[5] <- "TRT"
remove(heads)





##### Get mean and DS for initial rosette size distributions ----------------------------------------------

# Get distribution of CA rosette sizes; is normally distributed
fits_rs_CA <- fitdistr(na.omit(subset(data_rs, Species == "CA")$DM_t), "normal")$estimate
ks.test(na.omit(subset(data_rs, Species == "CA")$DM_t),
        pnorm, mean = fits_rs_CA[1], sd = fits_rs_CA[2])

# Get distribution of CN rosette sizes; is normally distributed
fits_rs_CN <- fitdistr(na.omit(subset(data_rs, Species == "CN")$DM_t), "normal")$estimate
ks.test(na.omit(subset(data_rs, Species == "CN")$DM_t),
        pnorm, mean = fits_rs_CN[1], sd = fits_rs_CN[2])





##### Get mean and SD for flower height distributions -----------------------------------------------------

# Create vector of CN and CA flower heights for each treatment group
ht_CN_NW <- subset(data_ht, Species == "CN" & TRT == "NW")
ht_CN_W <- subset(data_ht, Species == "CN" & TRT == "W")
ht_CA_NW <- subset(data_ht, Species == "CA" & TRT == "NW")
ht_CA_W <- subset(data_ht, Species == "CA" & TRT == "W")

# Get distribution of CA flower heights
fits_hd_CA_NW <- fitdistr(ht_CA_NW$Height, "normal")$estimate
ks.test(ht_CA_NW$Height, pnorm, mean = fits_hd_CA_NW[1], sd = fits_hd_CA_NW[2])
qqnorm(ht_CA_NW$Height)
qqline(ht_CA_NW$Height)
fits_hd_CA_W <- fitdistr(ht_CA_W$Height, "normal")$estimate
ks.test(ht_CA_W$Height, pnorm, mean = fits_hd_CA_W[1], sd = fits_hd_CA_W[2])
qqnorm(ht_CA_W$Height)
qqline(ht_CA_W$Height)

# Get distribution of CN flower heights
fits_hd_CN_NW <- fitdistr(ht_CN_NW$Height, "normal")$estimate
ks.test(ht_CN_NW$Height, pnorm, mean = fits_hd_CN_NW[1], sd = fits_hd_CN_NW[2])
qqnorm(ht_CN_NW$Height)
qqline(ht_CN_NW$Height)
fits_hd_CN_W <- fitdistr(ht_CN_W$Height, "normal")$estimate
ks.test(ht_CN_W$Height, pnorm, mean = fits_hd_CN_W[1], sd = fits_hd_CN_W[2])
qqnorm(ht_CN_W$Height)
qqline(ht_CN_W$Height)

# Remove variables that are no longer needed
remove(ht_CN_NW, ht_CN_W, ht_CA_NW, ht_CA_W)





##### Get equations for survival --------------------------------------------------------------------------

# Model survival rates as function of rosette size
glmer(Survival ~ DM_t + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Survival)))
glmer(Survival ~ DM_t + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Survival)))

# There are far too few deaths for logistic models to be reliable!
# Logistic regression or MLE would likely lead to small-sample bias
plot(subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Survival))$DM_T,
     subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Survival))$Survival)
plot(subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Survival))$DM_T,
     subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Survival))$Survival)

# Thus, rates will be simply be estimated as a constant
surv_rs_CA <- nrow(subset(data_rs, Species == "CA" & TRT == "NW" & Survival == 1))/
  nrow(subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Survival)))
surv_rs_CN <- nrow(subset(data_rs, Species == "CN" & TRT == "NW" & Survival == 1))/
  nrow(subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Survival)))





##### Get equations for flowering -------------------------------------------------------------------------

# Model flowering rates as function of rosette size
glmer(F ~ DM_t + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(F)))
glmer(F ~ DM_t + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(F)))

# There are far too few deaths for logistic models to be reliable!
# Logistic regression or MLE would likely lead to small-sample bias
plot(subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(F))$DM_t,
     subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(F))$F)
plot(subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(F))$DM_t,
     subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(F))$F)

# Thus, rates will be estimated independent of rosette size (i.e. as a constant)
flow_rs_CA <- nrow(subset(data_rs, Species == "CA" & TRT == "NW" & F == 1))/
  nrow(subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(F)))
flow_rs_CN <- nrow(subset(data_rs, Species == "CN" & TRT == "NW" & F == 1))/
  nrow(subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(F)))





##### Get equations for number of flower heads ------------------------------------------------------------

# Model CN number of flower heads as a function of rosette size
# Use AIC to make stepwise simplifications
mod_head_CN <- lmer(Heads ~ DM_t + (1|Row/Group),
                    data = subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Heads) & !is.na(DM_t)))
step(mod_head_CN)
mod_head_CN <- lmer(Heads ~ DM_t + (1|Group:Row),
                    data = subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Heads) & !is.na(DM_t)))
summary(mod_head_CN)
mod_head_CN <- fixef(mod_head_CN)

# Model CA number of flower heads as a function of rosette size
# Use AIC to make stepwise simplifications
mod_head_CA <- lmer(Heads ~ DM_t + (1|Row/Group),
                    data = subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Heads) & !is.na(DM_t)))
step(mod_head_CA)
mod_head_CA <- lmer(Heads ~ (1|Group:Row),
                    data = subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Heads) & !is.na(DM_t)))
summary(mod_head_CA)
mod_head_CA <- fixef(mod_head_CA)





##### Set up dispersal framework --------------------------------------------------------------------------

# Fit Weibull distribution to wind speeds
# Assume no seed release occurs for wind speeds of zero, so remove zero values
# Wind data is a bit messy, so KS is not significant; fit is still reasonable, though
ws_values <- c(data_ws1$Wind1, data_ws2$Wind1)
ws_values <- ws_values[ws_values > 0]
ws_params <- fitdistr(ws_values, "weibull")$estimate
ks.test(ws_values, pweibull, shape = ws_params[1], scale = ws_params[2])
plot(density(ws_values))
lines(density(rweibull(1000000, shape = ws_params[1], scale = ws_params[2])), col = "red")

# Function to transform raw variance and/or mean, then output corresponding shape and scale
transform.wb <- function(mean1, sd1, fv, which.trans){
  
  # Internal functions to calculate mean and SD of Weibull distribution
  meanW <- function(scale, shape){
    return(scale*gamma(1 + 1/shape))}
  sdW <- function(scale, shape){
    return(sqrt((scale^2)*(gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)))}
  
  # Must brute force a solution since there is no easy function for inverse gamma
  
  # Create meshpoints for combinations of shape and scale
  scale <- seq(0.1, 5, by = 0.005)
  shape <- seq(0.1, 5, by = 0.005)
  scale1 = rep(scale, each = length(shape))
  shape1 = rep(shape, times = length(scale))
  
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
  newMean = mean1*fm
  newSD = sd1*fs
  
  # Evaluate mean and SD for at meshpoints
  # Then calculate difference between estimated and supplied mean and SD
  diffs1 <- abs((mapply(meanW, scale = scale1, shape = shape1) - newMean)/newMean)
  diffs2 <- abs((mapply(sdW, scale = scale1, shape = shape1) - newSD)/newSD)
  
  # Find argmin by equally weighting mean and SD differences from supplied
  argmin <- which.min((diffs1 + diffs2)/2)
  
  # Return shape and scale at argmin
  return(c(shape1[argmin], scale1[argmin]))}

# Fit lognormal distribution to terminal velocities
# Terminal velocity is drop tube length (1.25 m) divided by drop time
tv_values_CN <- na.omit(1.25/subset(data_tv, species == "n")$drop.time)
tv_params_CN <- fitdistr(tv_values_CN, "lognormal")$estimate
ks.test(tv_values_CN, plnorm, meanlog = tv_params_CN[1], sdlog = tv_params_CN[2])
plot(density(tv_values_CN))
lines(density(rlnorm(1000000, meanlog = tv_params_CN[1], sdlog = tv_params_CN[2])), col = "red")

# Fit lognormal distribution to terminal velocities
# Terminal velocity is drop tube length (1.25 m) divided by drop time
tv_values_CA <- na.omit(1.25/subset(data_tv, species == "a")$drop.time)
tv_params_CA <- fitdistr(tv_values_CA, "lognormal")$estimate
ks.test(tv_values_CA, plnorm, meanlog = tv_params_CA[1], sdlog = tv_params_CA[2])
plot(density(tv_values_CA))
lines(density(rlnorm(1000000, meanlog = tv_params_CA[1], sdlog = tv_params_CA[2])), col = "red")

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

# Function generating a dispersal kernel using WALD model (Katul et al. 2005)
# Code adapted from Skarpaas and Shea (2007)
WALD.b <- function(n, H, species){
  
  # Prepare vector to store scalable dispersal parameters
  sParam <- c()
  
  # Set parameters for wind speed, seed terminal velocity, and vegetation height
  sParam[1] <- 0.15                 # Vegetation heignt in m
  sParam[2] <- ws_params[1]         # Mean wind speed, Weibull dist.
  sParam[3] <- ws_params[2]         # SD wind speed, Weibull dist.
  sParam[4] <- tv_params_CA[1]      # Mean terminal velocity, lognormal dist. (CA)
  sParam[5] <- tv_params_CA[2]      # SD terminal velocity, lognormal dist. (CA)
  sParam[6] <- tv_params_CN[1]      # Mean terminal velocity, lognormal dist. (CN)
  sParam[7] <- tv_params_CN[2]      # SD terminal velocity, lognormal dist. (CN)

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





##### Set up demography framework -------------------------------------------------------------------------

# Demography function for survival, reproduction, and more
demo <- function(dType, species, n = 0, rsize = 0, nflow = 0, dNum = dNum, dVal = dVal){
  
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
  dParam[14] <- NA                  # Prob. flowering slope, not applicable (CA)
  dParam[15] <- flow_rs_CN          # Prob. flowering intercept (CN)
  dParam[16] <- NA                  # Prob. flowering slope, not applicable (CN)
  dParam[17] <- mod_head_CA         # Num. heads intercept (CA)
  dParam[18] <- NA                  # Num. heads slope, not applicable (CA)
  dParam[19] <- mod_head_CN[1]      # Num. heads intercept (CN)
  dParam[20] <- mod_head_CN[2]      # Num. heads slope (CN)
  
  # Parameters for distribution of flower head heights among flowering individuals
  dParam[21] <- fits_hd_CA_NW[1]    # Mean head height (CA)
  dParam[22] <- fits_hd_CA_NW[2]    # SD head height (CA)
  dParam[23] <- fits_hd_CN_NW[1]    # Mean head height (CN)
  dParam[24] <- fits_hd_CN_NW[2]    # SD head height (CA)
  
  # Parameters for survival probability of rosettes
  dParam[25] <- surv_rs_CA          # Prob. survival intercept (CA)
  dParam[26] <- NA                  # Prob. survival slope, not applicable (CA)
  dParam[27] <- surv_rs_CN          # Prob. survival intercept (CN)
  dParam[28] <- NA                  # Prob. survival slope, not applicable (CN)
  
  # Scale only specified demographic parameter
  dVec <- rep(1, length(dParam))
  dVec[dNum] <- dVal
  dParam <- dParam*dVec
  
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
      prob <- dParam[13]}
    if(species == "CN"){
      prob <- dParam[15]}
    outcomes <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)
    return(outcomes)}
  
  # Flower production as function of initial rosette size
  # Round any non-integers up to the nearest head
  if(dType == "flowers"){
    if(species == "CA"){
      head <- rep(dParam[17], length(rsize))}
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
      prob <- dParam[25]}
    if(species == "CN"){
      prob <- dParam[27]}
    outcomes <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)}}





##### 1D expansion ----------------------------------------------------------------------------------------

# Function to simulate invasion wave
waveSim <- function(dNum, dVal){
  
  # Function to see if a seed is taken to the nearest nest
  nestsearch <- function(d, range){
    dists <- abs(d - nestsR)
    centre <- nestsR[which.min(dists)]
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
  nYear <- 1000     # Number of years to simulate
  trim <- TRUE      # Should core area of wave be trimmed?
  trimAmt <- 500    # Distance (m) behind wavefront to trim
  tDens <- 10       # Max thistle density per metre
  pAdult <- 0.5     # Proportion of rosettes that reach adulthood in 1 year
  plotOn <- FALSE   # Plot wave?
  
  # Run invasion wave simulation
  for(i in 1:nYear){
    
    # Initialise simulation data
    if(i == 1){
      
      # Generate nests
      nests <- sample(seq(0, nYear*200, by = 0.1), nDens*nYear*200) + 0.01
      
      # Initialise data; start with a single rosette
      plants <- data.frame(d = 0.01, stage = 0, rsize = demo("rsize", species, n = 1,
                                                             dNum = dNum, dVal = dVal), flow = 0)
      seedsAG <- data.frame(matrix(ncol = 2, nrow = 0))
      colnames(seedsAG) <- c("d", "germ")
      seedsSB <- data.frame(matrix(ncol = 2, nrow = 0))
      colnames(seedsSB) <- c("d", "germ")
      vals <- c()}
    
    # Surviving above-ground seeds become rosettes
    seedsAG$stage <- rep(0, nrow(seedsAG))
    seedsAG$rsize <- demo("rsize", species, n = nrow(seedsAG), dNum = dNum, dVal = dVal)
    seedsAG$flow <- rep(0, nrow(seedsAG))
    seedsAG <- seedsAG[, !names(seedsAG) == c("germ")]
    plants <- rbind(plants, seedsAG)
    
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
    
    # Reset seeds (except for seed bank)
    seedsAG <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsAG) <- c("d", "germ")
    seedsNew <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsNew) <- c("d", "germ")
    
    # Plot density over space
    if(i == 1){
      PlotList <- list()}
    if(plotOn == TRUE){
      PlotList[[i]] <- plants$d}
    
    # Remove rosettes that do not survive
    # Never kill when only one rosette; prevents code from crashing
    if(nrow(plants) > 1){
      plants <- plants[demo("survival", species, n = nrow(plants), dNum = dNum, dVal = dVal) == 1, ]}
    
    # Trim core areas as wave progresses to save computational resources
    if(trim == TRUE & nrow(plants) > 0){
      plants <- plants[plants$d > max(plants$d) - trimAmt, ]
      seedsSB <- seedsSB[seedsSB$d > max(plants$d) - trimAmt, ]
      nestsR <- nests[nests > min(plants$d) & nests < max(plants$d) + 1000]}
    
    # Some proportion of rosettes reach adulthood in one year
    # Ones that don't will remain rosettes for another year before becoming adults
    plants$stage[plants$stage == 1] <- 2
    plants$stage[plants$stage == 0] <- sample(c(1, 2), size = length(plants$stage[plants$stage == 0]),
                                              prob = c(1 - pAdult, pAdult), replace = TRUE)
    
    # Simulate flowering for adults
    if(sum(plants$stage == 2) > 0){
      plants$flow[plants$stage == 2] <- demo("flowering", species, dNum = dNum, dVal = dVal,
                                             n = length(plants$flow[plants$stage == 2]))}
    
    # Survival and establishment of seeds already in seed bank
    if(nrow(seedsSB) > 0){
      seedsSB <- seedsSB[demo("surSB", species, n = nrow(seedsSB), dNum = dNum, dVal = dVal) == 1, ]
      if(nrow(seedsSB) > 0){
        vec <- demo("estSB", species, n = nrow(seedsSB), dNum = dNum, dVal = dVal)
        temp1 <- seedsSB[vec == 1, ]
        temp1$stage <- rep(0, nrow(temp1))
        temp1$rsize <- demo("rsize", species, n = nrow(temp1), dNum = dNum, dVal = dVal)
        temp1$flow <- rep(0, nrow(temp1))
        temp1 <- temp1[, !names(temp1) == c("germ")]}
      plants <- na.omit(rbind(plants, temp1))
      seedsSB <- seedsSB[vec == 0, ]}
    
    # Simulate dispersal from flowering adults
    if(sum(plants$stage == 2) > 0){
      temp2 <- plants[plants$flow == 1, ]
      nflow <- demo("flowers", species, rsize = temp2$rsize, dNum = dNum, dVal = dVal)
      f1 <- unlist(sapply(nflow, demo, dType = "height", species = species, dNum = dNum, dVal = dVal))
      s1 <- rep(demo("seeds", "CN", dNum = dNum, dVal = dVal), times = sum(nflow))
      d1 <- unlist(mapply(x = temp2$d, times = nflow, rep))
      seedsNew <- as.vector(mapply(kern, n = s1, h = f1, d0 = d1, MoreArgs = list(species = species)))
      seedsNew <- data.frame(cbind(seedsNew, rep(0, length(seedsNew))))
      names(seedsNew) <- c("d", "germ")}
    
    # Simulate aboveground establishment from dispersed seeds
    seedsNew$germ <- demo("estAG", species, n = nrow(seedsNew), dNum = dNum, dVal = dVal)
    seedsAG <- na.omit(rbind(seedsAG, seedsNew[seedsNew$germ == 1, ]))
    seedsNew <- seedsNew[seedsNew$germ == 0, ]
    
    # Entry of non-establishing seeds into seed bank
    if(nrow(seedsNew) > 0){
      seedsNew$germ <- demo("entSB", species, n = nrow(seedsNew), dNum = dNum, dVal = dVal)
      seedsSB <- na.omit(rbind(seedsSB, seedsNew[seedsNew$germ == 1, ]))
      seedsSB$germ <- rep(0, nrow(seedsSB))}
    
    # Simulate secondary seed dispersal via ants
    if(nestOn == TRUE){
      seedsAG$d <- sapply(seedsAG$d, nestsearch, range = range)}
    
    # Store wavefront distance
    vals <- c(vals, max(plants$d))
    
    # Kill all adults after they reproduce
    plants <- plants[plants$stage != 2, ]}
  
  # Return list of wavefront positions and all plant positions
  return(list(wavefront = vals, positions = PlotList))}

# Calculate wavespeed elasticity
waveElas <- function(dNum, dVal){
  
  # Let dNum be the demographic parameter number
  # Let dVal be the proportion to multiply the original parameter by
  
  # Start timer
  time.start <- Sys.time()
  
  # Run simulation without any parameter changes
  wave1 <- waveSim(dNum = 1, dVal = 1)
  
  # Run simulation with increase/decrease on a specified parameter
  wave2 <- waveSim(dNum = dNum, dVal = dVal)
  
  # Calculate mean wavespeeds
  mWave1 <- mean(diff(wave1$wavefront))
  mWave2 <- mean(diff(wave2$wavefront))
  
  # Calculate percent change in wavespeed and parameter
  wPct <- (mWave2 - mWave1)/mWave1*100
  pPct <- (dVal - 1)*100
  
  # Calculate elasticity
  elas <- wPct/pPct
  
  # Calculate procedure time
  time.elapsed <- as.numeric(difftime(Sys.time(), time.start, units = "hours"))
  
  # Return list of calculated quantities
  return(list(procedure.time = paste0(time.elapsed, " hours"), elasticity = elas,
              wSpeedMean1 = mWave1, wSpeedMean2 = mWave2,
              wSpeed1 = diff(wave1$wavefront), wSpeed2 = diff(wave2$wavefront),
              wFront1 = wave1$wavefront, wFront2 = wave2$wavefront))}

# Calculate wavespeeds
#wave1 <- waveSim(dNum = 5, dVal = 0.7)
#mean(diff(wave1[[1]]))

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

