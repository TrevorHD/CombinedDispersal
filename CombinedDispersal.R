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
  
  # Simulate wind dispersal if released above canopy
  if(H > h){
    
    # Simulate wind speeds from empirical distribution of wind speeds
    Um <- rweibull(n, ws_params[1], ws_params[2])
  
    # Simulate terminal velocities from lognormal distribution
    if(species == "CN"){
      f <- rlnorm(n, tv_params_CN[1], tv_params_CN[2])}
    if(species == "CA"){
      f <- rlnorm(n, tv_params_CA[1], tv_params_CA[2])}
  
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
demo <- function(dType, species, n = 0, rsize = 0, nflow = 0){
  
  # Set various seed production and establishment parameters
  seedsCA <- 83       # Seeds per flower head (CA)
  seedsCN <- 160      # Seeds per flower head (CN)
  estCA <- 0.136      # Probability of establishment from seed (CA)
  estCN <- 0.147      # Probability of establishment from seed (CN)
  sPred <- 0.90       # Probability of seed predation
  sbEnt <- 0.233      # Probability of seed entering seed bank
  sbEst <- 0.233      # Probability of seed establishing from seed bank
  sbSur <- 0.260      # Probability of seed survival in seed bank
  
  # Per-head production of seeds, and subsequent seed survival
  if(dType == "seeds"){
    if(species == "CA"){
      nseed <- seedsCA*(1 - sPred)}
    if(species == "CN"){
      nseed <- seedsCN*(1 - sPred)}
    return(round(nseed))}
  
  # Establishment of seeds that do not enter the seed bank
  if(dType == "estAG"){
    if(species == "CA"){
      prob <- estCA}
    if(species == "CN"){
      prob <- estCN}
    outcomes <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)
    return(outcomes)}
  
  # Entry of seeds into the seed bank
  if(dType == "entSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(sbEnt, 1 - sbEnt), replace = TRUE)
    return(outcomes)}
  
  # Survival of seeds in the seed bank
  if(dType == "surSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(sbSur, 1 - sbSur), replace = TRUE)
    return(outcomes)}
  
  # Establishment of seeds from the seed bank
  if(dType == "estSB"){
    outcomes <- sample(c(1, 0), size = n, prob = c(sbEst, 1 - sbEst), replace = TRUE)
    return(outcomes)}
  
  # Initial rosette size upon establishment
  if(dType == "rsize"){
    if(species == "CA"){
      ros <- rtruncnorm(n, a = min(na.omit(data_rs$DM_t)),
                        mean = fits_rs_CA[1], sd = fits_rs_CA[2])}
    if(species == "CN"){
      ros <- rtruncnorm(n, a = min(na.omit(data_rs$DM_t)),
                        mean = fits_rs_CN[1], sd = fits_rs_CN[2])}
    return(ros)}
  
  # Flowering probability as function of initial rosette size
  if(dType == "flowering"){
    if(species == "CA"){
      prob <- flow_rs_CA}
    if(species == "CN"){
      prob <- flow_rs_CN}
    outcomes <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)
    return(outcomes)}
  
  # Flower production as function of initial rosette size
  # Round any non-integers up to the nearest head
  if(dType == "flowers"){
    if(species == "CA"){
      head <- rep(mod_head_CA, length(rsize))}
    if(species == "CN"){
      head <- mod_head_CN[1] + mod_head_CN[2]*rsize}
    return(ceiling(head))}
  
  # Distribution of flower heights for a given individual
  if(dType == "height"){
    if(species == "CA"){
      height <- rnorm(n, mean = fits_hd_CA_NW[1], sd = fits_hd_CA_NW[2])/100}
    if(species == "CN"){
      height <- rnorm(n, mean = fits_hd_CN_NW[1], sd = fits_hd_CN_NW[2])/100}
    return(height)}
  
  # Rosette survival as function of initial rosette size
  if(dType == "survival"){
    if(species == "CA"){
      prob <- surv_rs_CA}
    if(species == "CN"){
      prob <- surv_rs_CN}
    outcomes <- sample(c(1, 0), size = n, prob = c(prob, 1 - prob), replace = TRUE)}}





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
nYear <- 100      # Number of years to simulate
pAdult <- 0.5     # Proportion of rosettes that reach adulthood in 1 year
plotOn <- FALSE   # Plot wave?

# Run invasion wave simulation
for(i in 1:nYear){
  
  # Initialise simulation data
  if(i == 1){
    
    # Generate nests
    nests <- sample(seq(0, nYear*100, by = 0.1), nDens*nYear*100) + 0.01
    
    # Initialise data; start with a single rosette
    plants <- data.frame(d = 0.01, stage = 0, rsize = demo("rsize", species, n = 1), flow = 0)
    seedsAG <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsAG) <- c("d", "germ")
    seedsSB <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsSB) <- c("d", "germ")
    vals <- c()}
  
  # Surviving above-ground seeds become rosettes
  seedsAG$stage <- rep(0, nrow(seedsAG))
  seedsAG$rsize <- demo("rsize", species, n = nrow(seedsAG))
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
  if(plotOn == TRUE){
    if(i == 1){
      PlotList <- list()}
    PlotList[[i]] <- plants$d}
  
  # Remove rosettes that do not survive
  # Never kill when only one rosette; prevents code from crashing
  if(nrow(plants) > 1){
    plants <- plants[demo("survival", species, n = nrow(plants)) == 1, ]}
  
  # Trim core areas as wave progresses to save computational resources
  if(trim == TRUE & nrow(plants) > 0){
    plants <- plants[plants$d > max(plants$d) - trimAmt, ]
    seedsSB <- seedsSB[seedsSB$d > max(plants$d) - trimAmt, ]
    nests <- nests[nests > min(plants$d)]}
  
  # Some proportion of rosettes reach adulthood in one year
  # Ones that don't will remain rosettes for another year before becoming adults
  plants$stage[plants$stage == 1] <- 2
  plants$stage[plants$stage == 0] <- sample(c(1, 2), size = length(plants$stage[plants$stage == 0]),
                                            prob = c(1 - pAdult, pAdult), replace = TRUE)
  
  # Simulate flowering for adults
  if(sum(plants$stage == 2) > 0){
    plants$flow[plants$stage == 2] <- demo("flowering", species, n = length(plants$flow[plants$stage == 2]))}
  
  # Survival and establishment of seeds already in seed bank
  if(nrow(seedsSB) > 0){
    seedsSB <- seedsSB[demo("surSB", species, n = nrow(seedsSB)) == 1, ]
    if(nrow(seedsSB) > 0){
      vec <- demo("estSB", species, n = nrow(seedsSB))
      temp1 <- seedsSB[vec == 1, ]
      temp1$stage <- rep(0, nrow(temp1))
      temp1$rsize <- demo("rsize", species, n = nrow(temp1))
      temp1$flow <- rep(0, nrow(temp1))
      temp1 <- temp1[, !names(temp1) == c("germ")]}
    plants <- na.omit(rbind(plants, temp1))
    seedsSB <- seedsSB[vec == 0, ]}
  
  # Simulate dispersal from flowering adults
  if(sum(plants$stage == 2) > 0){
    temp2 <- plants[plants$flow == 1, ]
    nflow <- demo("flowers", species, rsize = temp2$rsize)
    f1 <- unlist(sapply(nflow, demo, dType = "height", species = species))
    s1 <- rep(demo("seeds", "CN"), times = sum(nflow))
    d1 <- unlist(mapply(x = temp2$d, times = nflow, rep))
    seedsNew <- as.vector(mapply(kern, n = s1, h = f1, d0 = d1, MoreArgs = list(species = species)))
    seedsNew <- data.frame(cbind(seedsNew, rep(0, length(seedsNew))))
    names(seedsNew) <- c("d", "germ")}
  
  # Simulate aboveground establishment from dispersed seeds
  seedsNew$germ <- demo("estAG", species, n = nrow(seedsNew))
  seedsAG <- na.omit(rbind(seedsAG, seedsNew[seedsNew$germ == 1, ]))
  seedsNew <- seedsNew[seedsNew$germ == 0, ]
  
  # Entry of non-establishing seeds into seed bank
  if(nrow(seedsNew) > 0){
    seedsNew$germ <- demo("entSB", species, n = nrow(seedsNew))
    seedsSB <- na.omit(rbind(seedsSB, seedsNew[seedsNew$germ == 1, ]))
    seedsSB$germ <- rep(0, nrow(seedsSB))}
  
  # Simulate secondary seed dispersal via ants
  if(nestOn == TRUE){
    seedsAG$d <- sapply(seedsAG$d, nestsearch, range = range)}
  
  # Store wavefront distance
  vals <- c(vals, max(plants$d))
  
  # Kill all adults after they reproduce
  plants <- plants[plants$stage != 2, ]}

# Get mean wavespeed
mean(diff(vals))

# Remove temporary variables
remove(temp1, temp2, vec, f1, s1, d1)

# Generate GIF of population spread
# Limit to 5000 m (nYear = 100 recommended)
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

