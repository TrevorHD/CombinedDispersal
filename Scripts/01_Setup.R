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

# Calculate plot averages for initial/final rosette size, survival, flowering, and head count
data_rs %>% 
  group_by(Row, Group, TRT) %>% 
  summarise(DM_t_PA = mean(DM_t, na.rm = TRUE),
            DM_t1_PA = mean(DM_t1, na.rm = TRUE),
            Survival_PA = mean(Survival, na.rm = TRUE),
            Flowering_PA = mean(Flowering, na.rm = TRUE),
            Heads_PA = mean(Heads, na.rm = TRUE)) -> data_rs_PA

# Create function to calculate rosette area
area <- function(diam){
  return(log(pi*(diam/2)^2))}





##### Get mean and SD for initial rosette size distributions ----------------------------------------------

# Model rosette area at t1 as a function of rosette area at t0; include warming and interaction
# Then perform stepwise selection to minimise AIC
mod_rose <- lmer(area(DM_t_PA) ~ TRT + (1|Group), data = data_rs_PA)
step(mod_rose)

# Dropping iwarming term minimises AIC; initial area is normally-distributed constant
mod_rose <- lmer(area(DM_t_PA) ~ (1|Group), data = data_rs_PA)
summary(mod_rose)

# Normally-distributed constant seems like a suitable approximation
plot(density(data_rs_PA$DM_t_PA))

# Model assumes that residuals are normally distributed around zero; assumption holds up
# Shapiro test can be sensitive to distribution tails, so take p-value witha grain of salt
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
ks.test(resid(mod_rose), pnorm, mean = 0, sd = sd(resid(mod_rose)))
shapiro.test(resid(mod_rose))
plot(density(resid(mod_rose)))
qqnorm(resid(mod_rose))
qqline(resid(mod_rose))
plot(mod_rose)

# Store model coefficients
# Get SD of errors to use as stochastic element in demographic simulations
mod_rose_err <- sd(resid(mod_rose))
mod_rose <- fixef(mod_rose)





##### Get equations for growth ----------------------------------------------------------------------------

# Model rosette area at t1 as a function of rosette area at t0; include warming and interaction
# Then perform stepwise selection to minimise AIC
mod_grow <- lmer(area(DM_t1_PA) ~ area(DM_t_PA) + TRT + TRT:area(DM_t_PA) + (1|Group), data = data_rs_PA)
step(mod_grow)

# Dropping interaction term minimises AIC
mod_grow <- lmer(area(DM_t1_PA) ~ area(DM_t_PA) + TRT + (1|Group), data = data_rs_PA)
summary(mod_grow)

# Model assumes that residuals are normally distributed around zero; assumption holds up
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
ks.test(resid(mod_grow), pnorm, mean = 0, sd = sd(resid(mod_grow)))
shapiro.test(resid(mod_grow))
plot(density(resid(mod_grow)))
qqnorm(resid(mod_grow))
qqline(resid(mod_grow))
plot(mod_grow)

# Store model coefficients
# Get SD of errors to use as stochastic element in demographic simulations
mod_grow_err <- sd(resid(mod_grow))
mod_grow_NW <- fixef(mod_grow)[c(1, 2)]
mod_grow_W <- fixef(mod_grow)[c(1, 2)] + c(fixef(mod_grow)[3], 0)
remove(mod_grow)





##### Get equations for survival --------------------------------------------------------------------------

# There are far too few deaths for logistic models to be reliable!
# Logistic regression or MLE would likely lead to small-sample bias
plot(subset(data_rs_PA, TRT == "NW")$DM_t_PA, subset(data_rs_PA, TRT == "NW")$Survival_PA)
plot(subset(data_rs_PA, TRT == "W")$DM_t_PA, subset(data_rs_PA, TRT == "W")$Survival_PA)

# Thus, rates will be simply be estimated as a constant
surv_rs_NW <- mean(subset(data_rs_PA, TRT == "NW")$Survival_PA)
surv_rs_W <- mean(subset(data_rs_PA, TRT == "W")$Survival_PA)





##### Get equations for flowering -------------------------------------------------------------------------

# There are far too few deaths for logistic models to be reliable!
# Logistic regression or MLE would likely lead to small-sample bias
plot(subset(data_rs_PA, TRT == "NW")$DM_t1_PA, subset(data_rs_PA, TRT == "NW")$Flowering_PA)
plot(subset(data_rs_PA, TRT == "W")$DM_t1_PA, subset(data_rs_PA, TRT == "W")$Flowering_PA)

# Thus, rates will be simply be estimated as a constant
flow_rs_NW <- mean(subset(data_rs_PA, TRT == "NW")$Flowering_PA)
flow_rs_W <- mean(subset(data_rs_PA, TRT == "W")$Flowering_PA)





##### Get equations for number of flower heads ------------------------------------------------------------

# Model flower head count as a function of rosette area at t0; include warming and interaction
# Then perform stepwise selection to minimise AIC
mod_head <- lmer(Heads_PA ~ area(DM_t1_PA) + TRT + TRT:area(DM_t1_PA) + (1|Group), data = data_rs_PA)
step(mod_head)

# Dropping interaction term minimises AIC
mod_head <- lmer(Heads_PA ~ area(DM_t1_PA) + TRT + (1|Group), data = data_rs_PA)
summary(mod_head)

# Model assumes that residuals are normally distributed around zero; assumption holds up
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
ks.test(resid(mod_head), pnorm, mean = 0, sd = sd(resid(mod_head)))
shapiro.test(resid(mod_head))
plot(density(resid(mod_head)))
qqnorm(resid(mod_head))
qqline(resid(mod_head))
plot(mod_head)

# Store model coefficients
# Get SD of errors to use as stochastic element in demographic simulations
mod_head_err <- sd(resid(mod_head))
mod_head_NW <- fixef(mod_head)[c(1, 2)]
mod_head_W <- fixef(mod_head)[c(1, 2)] + c(fixef(mod_head)[3], 0)
remove(mod_head)





##### Get mean and SD for flower height distributions -----------------------------------------------------

# Create vector of flower heights for each treatment group
ht_NW <- subset(data_ht, TRT == "NW")
ht_W <- subset(data_ht, TRT == "W")

# Get distribution of flower heights; is approximately normally distributed
fits_hd_NW <- fitdistr(ht_NW$Height, "normal")$estimate
ks.test(ht_NW$Height, pnorm, mean = fits_hd_NW[1], sd = fits_hd_NW[2])
shapiro.test(ht_NW$Height)
qqnorm(ht_NW$Height)
qqline(ht_NW$Height)
fits_hd_W <- fitdistr(ht_W$Height, "normal")$estimate
ks.test(ht_W$Height, pnorm, mean = fits_hd_W[1], sd = fits_hd_W[2])
shapiro.test(ht_W$Height)
qqnorm(ht_W$Height)
qqline(ht_W$Height)

# Get minimum and maximum observed height
ht_min <- min(data_ht$Height)
ht_max <- max(data_ht$Height)

# Remove variables that are no longer needed
remove(ht_NW, ht_W)





##### Set up wind speed and terminal velocity distributions -----------------------------------------------

# Fit Weibull distribution to wind speeds
# Assume no seed release occurs for wind speeds of zero, so remove zero values
# Wind data is a bit messy, so K-S is not significant; fit is still reasonable, though
ws_values <- c(data_ws1$Wind1, data_ws2$Wind1)
ws_values <- ws_values[ws_values > 0]
ws_params <- fitdistr(ws_values, "weibull")$estimate
ks.test(ws_values, pweibull, shape = ws_params[1], scale = ws_params[2])
plot(density(ws_values))
lines(density(rweibull(1000000, shape = ws_params[1], scale = ws_params[2])), col = "red")

# Fit lognormal distribution to terminal velocities
# Terminal velocity is drop tube length (1.25 m) divided by drop time
tv_values <- na.omit(1.25/subset(data_tv, species == "n")$drop.time)
tv_params <- fitdistr(tv_values, "lognormal")$estimate
ks.test(tv_values, plnorm, meanlog = tv_params[1], sdlog = tv_params[2])
plot(density(tv_values))
lines(density(rlnorm(1000000, meanlog = tv_params[1], sdlog = tv_params[2])), col = "red")





##### Set up ant seed dispersal distribution --------------------------------------------------------------

# Equation from Rabelo et al. (2021)... p = log(B0 + B1x)
# Probability of seed removal as function of distance from nest
an_params <- c(-0.56, -0.36)





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

