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
library(sqldf)
library(grid)
library(gridBase)
library(parallel)

# Set working directory
setwd("~/GitHub/CombinedDispersal")

# Load in raw weather data
data_ws1 <- read.csv("Data/Weather1.csv")
data_ws2 <- read.csv("Data/Weather2.csv")

# Load in terminal velocity data
# Use ambient TVs since we're only examining warming effects on height distribution
data_tv <- read.delim("Data/SeedDropData2.txt")
data_tv <- subset(data_tv, !is.na(drop.time))

# Load height data as xlsx
data_ht <- read.xlsx("Data/ThistleData.xlsx", sheetName = "Flowers")
data_ht <- subset(data_ht, (Type == "f" | Type == "s") & Species == "CN" & TRT != "PW",
                  select = -Type)

# Load rosette as xlsx
data_rs <- read.xlsx("Data/ThistleData.xlsx", sheetName = "General")
data_rs <- subset(data_rs, Species == "CN" & TRT != "PW",
                  select = -c(LLL_t, LLL_t1, OTC.On, OTC.Off, OTC.Notes))

# Rename rosette flowering column
names(data_rs)[9] <- "Flowering"

# Add survival indicator for rosettes that survived from establishment to summer
# Survival for rosettes that did not establish will be listed as NA
data_rs$Survival <- NA
data_rs$Survival[!is.na(data_rs$DM_t) & !is.na(data_rs$DM_t1)] <- 1
data_rs$Survival[!is.na(data_rs$DM_t) & is.na(data_rs$DM_t1)] <- 0

# Get number of flowers per plant
data_ht %>% 
  group_by(Row, Group, Plant, TRT) %>% 
  summarise(Heads = n()) %>% 
  data.frame() -> heads
data_rs <- merge(data_rs, heads, by = c("Row", "Group", "Plant"), all = TRUE)
data_rs <- subset(data_rs, select = -c(TRT.y))
names(data_rs)[5] <- "TRT"
remove(heads)

# Calculate plot averages for rosette size, survival, flowering, and head count
# Done to avoid pseudoreplication when conducting statistical analyses
data_rs %>% 
  group_by(Row, Group, TRT) %>% 
  summarise(DM_t_PA = mean(DM_t, na.rm = TRUE),
            DM_t1_PA = mean(DM_t1, na.rm = TRUE),
            Survival_PA = mean(Survival, na.rm = TRUE),
            Flowering_PA = mean(Flowering, na.rm = TRUE),
            Heads_PA = mean(Heads, na.rm = TRUE)) -> data_rs_PA

# Merge height data with rosette size data
data_ht <- merge(data_ht, data_rs, by = c("Row", "Group", "Plant"), all = TRUE)
data_ht <- subset(data_ht, select = c(Row, Group, Plant, Species.x, TRT.x, DM_t1, Height))
names(data_ht)[c(4, 5)] <- c("Species", "TRT")

# Calculate plot averages for flower height
# Done to avoid pseudoreplication when conducting statistical analyses
data_ht %>% 
  na.omit() %>% 
  group_by(Row, Group, Species, TRT) %>% 
  summarise(Height_PA = mean(Height)) -> data_ht_PA
data_ht_PA$DM_t1_PA <- data_rs_PA$DM_t1_PA

# Create function to calculate rosette area from diameter
# Also create inverse function that does the opposite
area <- function(diam){
  return(log(pi*(diam/2)^2))}
area.i <- function(area){
  return(2*sqrt(exp(area)/pi))}





##### Get mean and SD for initial rosette size distributions ----------------------------------------------

# Model initial rosette size at establishment as function of warming treatment
demo_rose <- lmer(area(DM_t_PA) ~ TRT + (1|Group), data = data_rs_PA)

# Then perform stepwise selection to minimise AIC
step(demo_rose)

# Dropping warming term minimises AIC
demo_rose <- lmer(area(DM_t_PA) ~ (1|Group), data = data_rs_PA)
summary(demo_rose)

# Model assumes residuals are normal around zero; assumption holds up
# Shapiro test is sensitive to distribution tails, so take p-value with grain of salt
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
ks.test(resid(demo_rose), pnorm, mean = 0, sd = sd(resid(demo_rose)))
shapiro.test(resid(demo_rose))
plot(density(resid(demo_rose)))
qqnorm(resid(demo_rose))
qqline(resid(demo_rose))
plot(demo_rose)

# Variances of error terms do not differ significantly between treatment groups
# We can thus model the error as normal with mean zero and SD agnostic of treatment
plot(density(resid(demo_rose)[which(data_rs_PA$TRT == "W")]))
lines(density(resid(demo_rose)[which(data_rs_PA$TRT == "NW")]), col = "red")
ks.test(resid(demo_rose)[which(data_rs_PA$TRT == "W")],
        resid(demo_rose)[which(data_rs_PA$TRT == "NW")])
var.test(resid(demo_rose)[which(data_rs_PA$TRT == "W")],
         resid(demo_rose)[which(data_rs_PA$TRT == "NW")])

# Store coefficients, and SD of errors to use as stochastic element in demographic simulations
demo_rose_err <- sd(resid(demo_rose))
demo_rose <- fixef(demo_rose)

# Test model fit; seems reasonable
set.seed(284759322)
temp1 <- area(data_rs_PA$DM_t_PA)
temp2 <- rep(demo_rose, length(temp1)) + rnorm(length(temp1), mean = 0, sd = demo_rose_err)
ks.test(temp1, temp2)
plot(density(temp1))
lines(density(temp2), col = "red")

# Remove unused variables
remove(temp1, temp2)





##### Get equations for survival --------------------------------------------------------------------------

# There are far too few deaths for logistic models to be reliable!
# Logistic regression or MLE would likely lead to small-sample bias
plot(subset(data_rs_PA, TRT == "NW")$DM_t_PA, subset(data_rs_PA, TRT == "NW")$Survival_PA)
plot(subset(data_rs_PA, TRT == "W")$DM_t_PA, subset(data_rs_PA, TRT == "W")$Survival_PA)

# Thus, rates will be simply be estimated as a single constant
demo_surv <- mean(data_rs_PA$Survival_PA)





##### Get equations for growth ----------------------------------------------------------------------------

# Model rosette area at t1 as a function of rosette area at t0; include warming and interaction
demo_grow <- lmer(area(DM_t1_PA) ~ area(DM_t_PA) + TRT + TRT:area(DM_t_PA) + (1|Group), data = data_rs_PA)

# Then perform stepwise selection to minimise AIC
step(demo_grow)

# Dropping interaction term minimises AIC
demo_grow <- lmer(area(DM_t1_PA) ~ area(DM_t_PA) + TRT + (1|Group), data = data_rs_PA)
summary(demo_grow)

# Model assumes residuals are normal around zero; assumption holds up
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
ks.test(resid(demo_grow), pnorm, mean = 0, sd = sd(resid(demo_grow)))
shapiro.test(resid(demo_grow))
plot(density(resid(demo_grow)))
qqnorm(resid(demo_grow))
qqline(resid(demo_grow))
plot(demo_grow)

# Store model coefficients
demo_grow_NW <- fixef(demo_grow)[c(1, 2)]
demo_grow_W <- fixef(demo_grow)[c(1, 2)] + c(fixef(demo_grow)[3], 0)

# Variances of error terms do not differ significantly between treatment groups
# We can thus model the error as normal with mean zero and SD agnostic of treatment
plot(density(resid(demo_grow)[which(data_rs_PA$TRT == "W")]))
lines(density(resid(demo_grow)[which(data_rs_PA$TRT == "NW")]), col = "red")
ks.test(resid(demo_grow)[which(data_rs_PA$TRT == "W")],
        resid(demo_grow)[which(data_rs_PA$TRT == "NW")])
var.test(resid(demo_grow)[which(data_rs_PA$TRT == "W")],
         resid(demo_grow)[which(data_rs_PA$TRT == "NW")])

# Store SD of errors to use as stochastic element in demographic simulations
demo_grow_err <- sd(resid(demo_grow))

# Test model fit (unwarmed); seems reasonable
set.seed(386589364)
temp1 <- area(subset(data_rs_PA, TRT == "NW")$DM_t1_PA)
temp2 <- demo_grow_NW[1] + demo_grow_NW[2]*area(subset(data_rs_PA, TRT == "NW")$DM_t_PA) +
  rnorm(length(temp1), mean = 0, sd = demo_grow_err)
ks.test(temp1, temp2)
plot(density(temp1))
lines(density(temp2), col = "red")

# Test model fit (warmed); seems reasonable
set.seed(386589364)
temp1 <- area(subset(data_rs_PA, TRT == "W")$DM_t1_PA)
temp2 <- demo_grow_W[1] + demo_grow_W[2]*area(subset(data_rs_PA, TRT == "W")$DM_t_PA) +
  rnorm(length(temp1), mean = 0, sd = demo_grow_err)
ks.test(temp1, temp2)
plot(density(temp1))
lines(density(temp2), col = "red")

# Remove unused variables
remove(demo_grow, temp1, temp2)





##### Get equations for flowering -------------------------------------------------------------------------

# There are far too few deaths for logistic models to be reliable!
# Logistic regression or MLE would likely lead to small-sample bias
plot(subset(data_rs_PA, TRT == "NW")$DM_t1_PA, subset(data_rs_PA, TRT == "NW")$Flowering_PA)
plot(subset(data_rs_PA, TRT == "W")$DM_t1_PA, subset(data_rs_PA, TRT == "W")$Flowering_PA)

# Thus, rates will be simply be estimated as a single constant
demo_flow <- mean(data_rs_PA$Flowering_PA)





##### Get equations for number of flower heads ------------------------------------------------------------

# Model flower head count as a function of rosette area at t1; include warming and interaction
demo_head <- lmer(Heads_PA ~ area(DM_t1_PA) + TRT + TRT:area(DM_t1_PA) + (1|Group), data = data_rs_PA)

# Perform stepwise selection to minimise AIC
step(demo_head)

# Dropping interaction term minimises AIC
demo_head <- lmer(Heads_PA ~ area(DM_t1_PA) + TRT + (1|Group), data = data_rs_PA)
summary(demo_head)

# Model assumes residuals are normal around zero; assumption holds up
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
ks.test(resid(demo_head), pnorm, mean = 0, sd = sd(resid(demo_head)))
shapiro.test(resid(demo_head))
plot(density(resid(demo_head)))
qqnorm(resid(demo_head))
qqline(resid(demo_head))
plot(demo_head)

# Store model coefficients
demo_head_NW <- fixef(demo_head)[c(1, 2)]
demo_head_W <- fixef(demo_head)[c(1, 2)] + c(fixef(demo_head)[3], 0)

# Variances of error terms do not differ significantly between treatment groups
# We can thus model the error as normal with mean zero and SD agnostic of treatment
plot(density(resid(demo_head)[which(data_rs_PA$TRT == "W")]))
lines(density(resid(demo_head)[which(data_rs_PA$TRT == "NW")]), col = "red")
ks.test(resid(demo_head)[which(data_rs_PA$TRT == "W")],
        resid(demo_head)[which(data_rs_PA$TRT == "NW")])
var.test(resid(demo_head)[which(data_rs_PA$TRT == "W")],
         resid(demo_head)[which(data_rs_PA$TRT == "NW")])

# Store SD of errors to use as stochastic element in demographic simulations
demo_head_err <- sd(resid(demo_head))

# Test model fit (unwarmed); seems reasonable
set.seed(927494737)
temp1 <- subset(data_rs_PA, TRT == "NW")$Heads_PA
temp2 <- demo_head_NW[1] + demo_head_NW[2]*area(subset(data_rs_PA, TRT == "NW")$DM_t1_PA) +
  rnorm(length(temp1), mean = 0, sd = demo_head_err)
ks.test(temp1, temp2)
plot(density(temp1))
lines(density(temp2), col = "red")

# Test model fit (warmed); seems reasonable
set.seed(927494737)
temp1 <- subset(data_rs_PA, TRT == "W")$Heads_PA
temp2 <- demo_head_W[1] + demo_head_W[2]*area(subset(data_rs_PA, TRT == "W")$DM_t1_PA) +
  rnorm(length(temp1), mean = 0, sd = demo_head_err)
ks.test(temp1, temp2)
plot(density(temp1))
lines(density(temp2), col = "red")

# Remove unused variables
remove(demo_head, temp1, temp2)





##### Get equations for flower head heights ---------------------------------------------------------------

# Model flower head height as a function of rosette size; include warming and interaction
# Keep structure consistent with Drees and Shea (2023) by using diameter as covariate instead of area
# However, use diameter at t1 instead of t0; makes more sense in our Demo + Dispersal model framework
# Thus, parameters estimates will be slightly different compared to Drees and Shea (2023)
demo_hdht <- lmer(Height_PA ~ DM_t1_PA + TRT + TRT:DM_t1_PA + (1|Group), data = data_ht_PA)

# Perform stepwise selection to minimise AIC
step(demo_hdht)

# Dropping interaction term minimises AIC
demo_hdht <- lmer(Height_PA ~ DM_t1_PA + TRT + (1|Group), data = data_ht_PA)
summary(demo_hdht)

# Model assumes residuals are normal around zero; assumption holds up
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
ks.test(resid(demo_hdht), pnorm, mean = 0, sd = sd(resid(demo_hdht)))
shapiro.test(resid(demo_hdht))
plot(density(resid(demo_hdht)))
qqnorm(resid(demo_hdht))
qqline(resid(demo_hdht))
plot(demo_hdht)

# Store model coefficients
demo_hdht_NW <- fixef(demo_hdht)[c(1, 2)]
demo_hdht_W <- fixef(demo_hdht)[c(1, 2)] + c(fixef(demo_hdht)[3], 0)

# Estimate errors, but use non-PA distribution since plot averaging mutes small and large heights
# This drastically shrinks the variance of the true distribution of flower head heights
# In turn, this could significantly over- or under-estimate spread rate
temp1 <- drop_na(data_ht)
names(temp1)[6] <- c("DM_t1_PA")
temp2 <- predict(demo_hdht, temp1) - temp1$Height

# Variances of error terms do not differ significantly between treatment groups
# We can thus model the error as normal with mean zero and SD agnostic of treatment
plot(density(temp2[which(temp1$TRT == "W")]))
lines(density(temp2[which(temp1$TRT == "NW")]), col = "red")
ks.test(temp2[which(temp1$TRT == "W")], temp2[which(temp1$TRT == "NW")])
var.test(temp2[which(temp1$TRT == "W")], temp2[which(temp1$TRT == "NW")])

# Store SD of errors to use as stochastic element in demographic simulations
demo_hdht_err <- sd(temp2)

# Test model fit (unwarmed); seems reasonable
set.seed(284657777)
temp3 <- drop_na(subset(data_ht, TRT == "NW"))$Height
temp4 <- demo_hdht_NW[1] + demo_hdht_NW[2]*drop_na(subset(data_ht, TRT == "NW"))$DM_t1 +
  rnorm(length(temp3), mean = 0, sd = demo_hdht_err)
ks.test(temp3, temp4)
plot(density(temp3))
lines(density(temp4), col = "red")

# Test model fit (warmed); seems reasonable
set.seed(284657777)
temp5 <- drop_na(subset(data_ht, TRT == "W"))$Height
temp6 <- demo_hdht_W[1] + demo_hdht_W[2]*drop_na(subset(data_ht, TRT == "W"))$DM_t1 +
  rnorm(length(temp5), mean = 0, sd = demo_hdht_err)
ks.test(temp5, temp6)
plot(density(temp5))
lines(density(temp6), col = "red")

# Remove unused variables
remove(demo_hdht, temp1, temp2, temp3, temp4, temp5, temp6)





##### Set up wind speed and terminal velocity distributions -----------------------------------------------

# Fit Weibull distribution to wind speeds
# Assume no seed release occurs for wind speeds of zero, so remove zero values
disp_ws_vals <- c(data_ws1$Wind1, data_ws2$Wind1)
disp_ws_vals <- disp_ws_vals[disp_ws_vals > 0]
disp_ws <- fitdistr(disp_ws_vals, "weibull")$estimate

# Get mean and SD for parameterised Weibull distribution
disp_ws[2]*gamma(1 + 1/disp_ws[1])
sqrt((disp_ws[2]^2)*(gamma(1 + 2/disp_ws[1]) - gamma(1 + 1/disp_ws[1])^2))

# Wind data is a bit messy, so K-S is not significant; fit is still reasonable, though
# Note that previous C. nutans studies also use Weibull for wind speeds
set.seed(283749842)
temp1 <- rweibull(length(disp_ws_vals), shape = disp_ws[1], scale = disp_ws[2])
ks.test(disp_ws_vals, temp1)
plot(density(disp_ws_vals))
lines(density(temp1), col = "red")

# Fit lognormal distribution to terminal velocities
# Terminal velocity is drop tube length (1.25 m) divided by drop time
disp_tv_vals <- na.omit(1.25/subset(data_tv, species == "n")$drop.time)
disp_tv <- fitdistr(disp_tv_vals, "lognormal")$estimate

# Get mean and SD for parameterised lognormal distribution
exp(disp_tv[1] + 0.5*(disp_tv[2]^2))
sqrt(exp(2*disp_tv[1] + (disp_tv[2]^2))*(exp(disp_tv[2]^2) - 1))

# Lognormal seems to be reasonable fit for terminal velocity data
set.seed(283749842)
temp2 <- rlnorm(length(disp_tv_vals), meanlog = disp_tv[1], sdlog = disp_tv[2])
ks.test(disp_tv_vals, temp2)
plot(density(disp_tv_vals))
lines(density(temp2), col = "red")

# Remove unused variables
remove(temp1, temp2)





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





##### Set up functions for density calculations -----------------------------------------------------------

# Function to calculate density per 1-m window
dens.calc <- function(i, plants){
  
  # Get list of positions within a 1-m window
  vals <- plants$d[plants$d >= i - 1 & plants$d < i]
  
  # Get number of plants in the same window
  dens <- length(vals)
  
  # Return density
  return(dens)}

# Function to spread plants evenly in max-density windows
# More realistic than leaving all plants in exact same position
dens.spread <- function(i, plants){
  
  # Get list of plants within a 1-m window
  vals <- plants[plants$d >= i - 1 & plants$d < i, ]
  
  # Spread plants evenly if at max density
  if(nrow(vals) == tDens){
    vals$d <- i - 1 + (0:(tDens - 1))/tDens}
  
  # Return updated dataframe of plants
  return(vals)}





##### Create function to plot wavefronts over time --------------------------------------------------------

# Function to generate GIF of population spread (nYear = 100 recommended)
generatePlots <- function(type = "hist", snap = FALSE, snap_t = NULL, snap_bottom = FALSE){
  
  # Set x-axis maximum based on wavefront at end of simulation (or a specified number if using snapshots)
  # Set other graphical parameters based on whether or not there will be multiple adjacent plots
  if(snap == FALSE){
    max_scale <- (ceiling(max(wv_plots[[length(wv_plots)]])/1000) + 1)*1000}
  if(snap == TRUE){
    max_scale <- 3000}
  xlab <- ifelse(snap_bottom == TRUE | snap == FALSE, "Distance (m)", "")
  tsize <- ifelse(snap == FALSE, 1, 0.5)
  
  if(snap == FALSE){
    plotSeq <- 1:length(wv_plots)}
  if(snap == TRUE){
    plotSeq <- snap_t}
  
  # Generate plot using position data at each timestep
  for(i in plotSeq){
    
    # Initialise vectors
    lower <- c()
    positions <- wv_plots[[i]]
    if(length(positions) == 0){
      positions <- 0}
    
    # Assume max density behind trimmed wavefront
    # Trim condition here should match that used in the simulation
    if(max(positions > 1000)){
      minBin <- sum(positions > floor(min(positions)) & positions < ceiling(min(positions)))
      lower <- c(rep(0:(floor(min(positions) - 1)), 10) + 0.01,
                 rep(floor(min(positions)) + 0.01, 10 - minBin))
      if(min(positions) < 1){
        lower <- sort(lower)[-c(1:20)]}}
    
    # Plot wave movement as either histogram or density relative to maximum
    if(type == "hist"){
      hist(c(lower, positions), breaks = seq(0, max_scale, by = 20), xlim = c(0, max_scale), ylim = c(0, 215),
           xaxt = "n", yaxt = "n", xlab = NULL, ylab = "Count", main = "")
      axis(2, at = seq(0, 200, 50))}
    if(type == "density"){
      plotdata <- hist(c(lower, positions), breaks = seq(0, max_scale, by = 20), plot = FALSE)
      plot(plotdata$mids, plotdata$density/max(plotdata$density),
           type = "l", xlim = c(0, max_scale), ylim = c(0, 1.075),
           xaxt = "n", yaxt = "n", xlab = "", ylab = "Relative Density", main = "")
      axis(2, at = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"))}
    if(snap == FALSE){
      axis(1, at = seq(0, max_scale, 1000))
      title(xlab = xlab)}
    if(snap == TRUE){
      if(snap_bottom == FALSE){
        axis(1, at = seq(0, max_scale, 1000), labels = FALSE)}
      if(snap_bottom == TRUE){
        axis(1, at = seq(0, max_scale, 1000), mgp = c(0, -0.1, 0))
        title(xlab = xlab, mgp = c(0.5, 0, 0))}}
    if(type == "hist"){
      text(x = max_scale*0.95, y = 200, paste0("t = ", i), cex = tsize)
      box()}
    if(type == "density"){
      text(x = max_scale*0.95, y = 1, paste0("t = ", i), cex = tsize)}}}

