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
data_rs <- subset(data_rs, Species == "CN", select = -c(TRT.y))
names(data_rs)[5] <- "TRT"
remove(heads)





##### Get mean and SD for initial rosette size distributions ----------------------------------------------

# Get distribution of rosette sizes; is normally distributed
fits_rs <- fitdistr(na.omit(data_rs$DM_t), "normal")$estimate
ks.test(na.omit(data_rs$DM_t), pnorm, mean = fits_rs[1], sd = fits_rs[2])





##### Get mean and SD for flower height distributions -----------------------------------------------------

# Create vector of flower heights for each treatment group
ht_NW <- subset(data_ht, TRT == "NW")
ht_W <- subset(data_ht, TRT == "W")

# Get distribution of flower heights
fits_hd_NW <- fitdistr(ht_NW$Height, "normal")$estimate
ks.test(ht_NW$Height, pnorm, mean = fits_hd_NW[1], sd = fits_hd_NW[2])
qqnorm(ht_NW$Height)
qqline(ht_NW$Height)
fits_hd_W <- fitdistr(ht_W$Height, "normal")$estimate
ks.test(ht_W$Height, pnorm, mean = fits_hd_W[1], sd = fits_hd_W[2])
qqnorm(ht_W$Height)
qqline(ht_W$Height)

# Get minimum and maximum observed height
ht_min <- min(subset(data_ht, Species == "CN", select = Height))
ht_max <- max(subset(data_ht, Species == "CN", select = Height))

# Remove variables that are no longer needed
remove(ht_NW, ht_W)





##### Get equations for survival --------------------------------------------------------------------------

# Model survival rates as function of rosette size (area)
glmer(Survival ~ log(pi*(DM_t/2)^2) + (1|Row/Group), family = "binomial",
      data = subset(data_rs, TRT == "NW" & !is.na(Survival)))

# There are far too few deaths for logistic models to be reliable!
# Logistic regression or MLE would likely lead to small-sample bias
plot(subset(data_rs, TRT == "NW" & !is.na(Survival))$DM_t,
     subset(data_rs, TRT == "NW" & !is.na(Survival))$Survival)

# Thus, rates will be simply be estimated as a constant
surv_rs <- nrow(subset(data_rs, TRT == "NW" & Survival == 1))/
           nrow(subset(data_rs, TRT == "NW" & !is.na(Survival)))





##### Get equations for flowering -------------------------------------------------------------------------

# Model flowering rates as function of rosette size (area)
glmer(F ~ log(pi*(DM_t/2)^2) + (1|Row/Group), family = "binomial",
      data = subset(data_rs, TRT == "NW" & !is.na(F)))

# There are far too few deaths for logistic models to be reliable!
# Logistic regression or MLE would likely lead to small-sample bias
plot(subset(data_rs, TRT == "NW" & !is.na(F))$DM_t,
     subset(data_rs, TRT == "NW" & !is.na(F))$F)

# Thus, rates will be estimated independent of rosette size (i.e. as a constant)
flow_rs <- nrow(subset(data_rs, TRT == "NW" & F == 1))/
           nrow(subset(data_rs, TRT == "NW" & !is.na(F)))





##### Get equations for number of flower heads ------------------------------------------------------------

# Model number of flower heads as a function of rosette size (area)
# Use AIC to make stepwise simplifications
mod_head <- lmer(Heads ~ log(pi*(DM_t/2)^2) + (1|Row/Group),
                 data = subset(data_rs, TRT == "NW" & !is.na(Heads) & !is.na(DM_t)))
step(mod_head)
summary(mod_head)
mod_head <- fixef(mod_head)





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

