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
data_rs <- subset(data_rs, !is.na(Species), select = -c(TRT.y))
names(data_rs)[5] <- "TRT"
remove(heads)





##### Get mean and SD for initial rosette size distributions ----------------------------------------------

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

# Model survival rates as function of rosette size (area)
glmer(Survival ~ log(pi*(DM_t/2)^2) + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Survival)))
glmer(Survival ~ log(pi*(DM_t/2)^2) + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Survival)))

# There are far too few deaths for logistic models to be reliable!
# Logistic regression or MLE would likely lead to small-sample bias
plot(subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Survival))$DM_t,
     subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Survival))$Survival)
plot(subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Survival))$DM_t,
     subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Survival))$Survival)

# Thus, rates will be simply be estimated as a constant
surv_rs_CA <- nrow(subset(data_rs, Species == "CA" & TRT == "NW" & Survival == 1))/
  nrow(subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Survival)))
surv_rs_CN <- nrow(subset(data_rs, Species == "CN" & TRT == "NW" & Survival == 1))/
  nrow(subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Survival)))





##### Get equations for flowering -------------------------------------------------------------------------

# Model flowering rates as function of rosette size (area)
glmer(F ~ log(pi*(DM_t/2)^2) + (1|Row/Group), family = "binomial",
      data = subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(F)))
glmer(F ~ log(pi*(DM_t/2)^2) + (1|Row/Group), family = "binomial",
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

# Model CN number of flower heads as a function of rosette size (area)
# Use AIC to make stepwise simplifications
mod_head_CN <- lmer(Heads ~ log(pi*(DM_t/2)^2) + (1|Row/Group),
                    data = subset(data_rs, Species == "CN" & TRT == "NW" & !is.na(Heads) & !is.na(DM_t)))
step(mod_head_CN)
summary(mod_head_CN)
mod_head_CN <- fixef(mod_head_CN)

# Model CA number of flower heads as a function of rosette size (area)
# Use AIC to make stepwise simplifications
mod_head_CA <- lmer(Heads ~ log(pi*(DM_t/2)^2) + (1|Row/Group),
                    data = subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Heads) & !is.na(DM_t)))
step(mod_head_CA)
mod_head_CA <- lmer(Heads ~ (1|Group:Row),
                    data = subset(data_rs, Species == "CA" & TRT == "NW" & !is.na(Heads) & !is.na(DM_t)))
summary(mod_head_CA)
mod_head_CA <- fixef(mod_head_CA)





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





##### Set up ant seed dispersal distribution --------------------------------------------------------------

# Equation from Rabelo et al. (2021)... p = log(B0 + B1x)
# Probability of seed removal as function of distance from nest
an_params <- c(-0.56, -0.36)

