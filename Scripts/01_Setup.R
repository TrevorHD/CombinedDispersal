##### Load packages and initialise data -------------------------------------------------------------------

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
data_tv <- read.delim("Data/SeedDropData2.txt")
data_tv <- subset(data_tv, !is.na(drop.time) & species == "n")

# Load height data
data_ht <- read.xlsx("Data/ThistleData.xlsx", sheetName = "Flowers")
data_ht <- subset(data_ht, (Type == "f" | Type == "s") & Species == "CN" & TRT != "PW",
                  select = -Type)

# Load rosette data
data_rs <- read.xlsx("Data/ThistleData.xlsx", sheetName = "General")
data_rs <- subset(data_rs, Species == "CN" & TRT != "PW",
                  select = -c(LLL_t, LLL_t1, OTC.On, OTC.Off, OTC.Notes))

# Rename bolting/flowering indicator
names(data_rs)[9] <- "Flowering"

# Add survival indicator for rosettes that survived over winter
# Survival for rosettes that did not establish will be listed as NA
data_rs$Survival <- NA
data_rs$Survival[!is.na(data_rs$DM_t) & !is.na(data_rs$DM_t1)] <- 1
data_rs$Survival[!is.na(data_rs$DM_t) & is.na(data_rs$DM_t1)] <- 0

# Count flower heads per plant
data_ht %>% 
  group_by(Row, Group, Plant, TRT) %>% 
  summarise(Heads = n()) %>% 
  data.frame() -> heads
data_rs <- merge(data_rs, heads, by = c("Row", "Group", "Plant"), all = TRUE)
data_rs <- subset(data_rs, select = -c(TRT.y))
names(data_rs)[5] <- "TRT"
remove(heads)

# Calculate plot averages for rosette size, survival, flowering/bolting, and head count
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

# Calculate plot averages for flower head height
# Done to avoid pseudoreplication when conducting statistical analyses
data_ht %>% 
  na.omit() %>% 
  group_by(Row, Group, Species, TRT) %>% 
  summarise(Height_PA = mean(Height)) -> data_ht_PA
data_ht_PA$DM_t1_PA <- data_rs_PA$DM_t1_PA





##### Set up terminal velocity and wind speed transformation functions ------------------------------------

# Define function to transform raw variance and/or mean, then output corresponding meanlog and sdlog
transform.ln <- function(meanlog, sdlog, fv, which.trans){
  
  # Define internal functions to convert between mean/meanlog and sd/sdlog
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
  # Note: fv cannot equal zero when adjusting the mean
  
  # Set conditional statements for the variable(s) to be transformed
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

# Define function to transform raw variance and/or mean, then output corresponding shape and scale
transform.wb <- function(shape, scale, fv, which.trans){
  
  # Define internal functions to calculate mean and SD of Weibull distribution
  meanW <- function(shape, scale){
    return(scale*gamma(1 + 1/shape))}
  sdW <- function(shape, scale){
    return(sqrt((scale^2)*(gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)))}
  
  # Unfortunately, inverse gamma function cannot be expressed using elementary functions
  # Thus, it might be easier just to brute force a solution using a gridsearch 
  
  # Create meshpoints for combinations of shape and scale
  mesh <- seq(0, 5, by = 0.002)
  if(fv >= 1){mesh2 <- mesh[mesh >= scale*0.8]}
  if(fv < 1){mesh2 <- mesh[mesh < scale*1.2]}
  shapeMesh = rep(mesh, times = length(mesh2))
  scaleMesh = rep(mesh2, each = length(mesh))
  
  # Let fv be the value that the mean and/or standard deviation is multiplied by
  # Note: transformations may not be accurate when fv < 0.01 or fv > 4
  
  # Set conditional statements for the variable(s) to be transformed
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
  
  # Evaluate mean and SD at meshpoints
  # Then calculate difference between estimated and supplied mean and SD
  diffs1 <- abs((mapply(meanW, shape = shapeMesh, scale = scaleMesh) - newMean)/newMean)
  diffs2 <- abs((mapply(sdW, shape = shapeMesh, scale = scaleMesh) - newSD)/newSD)
  
  # Find argmin by equally weighting mean and SD differences from supplied
  argmin <- which.min((diffs1 + diffs2)/2)
  
  # Return shape and scale at argmin
  return(c(shapeMesh[argmin], scaleMesh[argmin]))}





##### Set up functions for plotting -----------------------------------------------------------------------

# Define function to plot wavespeed values across a span of parameter values
plot.span <- function(x, y, left = FALSE, bottom = FALSE, type = "real"){
  
  # Set baseline graphical parameters, axis labels, and colours
  par(mar = c(0.85, 0.85, 0, 0), mgp = c(0.35, 0, 0), cex.axis = 0.3, cex.lab = 0.4, tcl = -0.12)
  if(type == "real"){
    xlab <- "Prop. of Base Value"
    colour <- "green4"}
  if(type == "dens"){
    xlab <- "Prop. of Base Value"
    colour <- "purple"}
  if(type == "prob"){
    xlab <- "Probability"
    colour <- "green"}
  
  # Plot wavespeeds spanning across several parameter values
  if(left == FALSE & bottom == FALSE){
    plot(x, y, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, 100),
         cex = 0.5, col = colour, pch = 16)
    axis(side = 1, at = x, labels = FALSE)
    axis(side = 2, at = seq(0, 100, by = 10), labels = FALSE)}
  if(left == TRUE & bottom == FALSE){
    plot(x, y, xaxt = "n", yaxt = "n", xlab = "", ylab = "Spread Rate (m/yr)", ylim = c(0, 100),
         cex = 0.5, col = colour, pch = 16)
    axis(side = 1, at = x, labels = FALSE)
    axis(side = 2, at = seq(0, 100, by = 10))}
  if(left == FALSE & bottom == TRUE){
    plot(x, y, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, 100),
         cex = 0.5, col = colour, pch = 16)
    axis(side = 1, at = x, mgp = c(0, -0.4, 0))
    mtext(side = 1, line = -0.03, xlab, cex = 0.4)
    axis(side = 2, at = seq(0, 100, by = 10), labels = FALSE)}
  if(left == TRUE & bottom == TRUE){
    plot(x, y, xaxt = "n", yaxt = "n", xlab = "", ylab = "Spread Rate (m/yr)", ylim = c(0, 100),
         cex = 0.5, col = colour, pch = 16)
    axis(side = 1, at = x, mgp = c(0, -0.4, 0))
    mtext(side = 1, line = -0.03, xlab, cex = 0.4)
    axis(side = 2, at = seq(0, 100, by = 10))}}

# Define function to plot population spread (sets_years=100, misc_tDens=10 required)
plot.wave <- function(type = "hist", snap = FALSE, snap_t = NULL, snap_bottom = FALSE){
  
  # Set x-axis maximum based on wavefront at end of simulation (or a specified number if using snapshots)
  # Set other graphical parameters based on whether or not there will be multiple adjacent plots
  if(snap == FALSE){
    max_scale <- (ceiling(max(wv_plots[[length(wv_plots)]])/1000) + 1)*1000}
  if(snap == TRUE){
    max_scale <- 1600}
  xlab <- ifelse(snap_bottom == TRUE | snap == FALSE, "Distance (m)", "")
  tsize <- ifelse(snap == FALSE, 1, 0.5)
  if(snap == FALSE){
    plotSeq <- 1:(length(wv_plots) - 1)}
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
    plotdata <- hist(c(lower, positions), breaks = seq(0, max_scale, by = 20), plot = FALSE)
    plotdata$counts[plotdata$counts > 200] <- 200
    if(type == "hist"){
      plot(plotdata, xlim = c(0, max_scale), ylim = c(0, 215), xaxt = "n", yaxt = "n",
           xlab = NULL, ylab = "Relative Density", main = "")
      axis(2, at = seq(0, 200, 50), labels = paste0(seq(0, 1, 0.25)*100, "%"))}
    if(type == "line"){
      plot(plotdata$mids, plotdata$counts/200,
           type = "l", xlim = c(0, max_scale), ylim = c(0, 1.075), xaxt = "n", yaxt = "n",
           xlab = "", ylab = "Relative Density", main = "")
      axis(2, at = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"))}
    if(snap == FALSE){
      axis(1, at = seq(0, max_scale, 1000))
      title(xlab = xlab)}
    if(snap == TRUE){
      if(snap_bottom == FALSE){
        axis(1, at = seq(0, max_scale, 200), labels = FALSE)}
      if(snap_bottom == TRUE){
        axis(1, at = seq(0, max_scale, 200), mgp = c(0, -0.1, 0))
        title(xlab = xlab, mgp = c(0.5, 0, 0))}}
    if(type == "hist"){
      text(x = max_scale*0.95, y = 200, paste0("t = ", i), cex = tsize)
      box()}
    if(type == "line"){
      text(x = max_scale*0.95, y = 1, paste0("t = ", i), cex = tsize)}}}

