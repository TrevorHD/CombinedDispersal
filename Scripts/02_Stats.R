##### Get mean and SD for initial rosette size distributions ----------------------------------------------

# Plot data; mean and SD seem similar between treatments
plot(density(data_rs_PA$DM_t_PA[which(data_rs_PA$TRT == "W")]), col = "red",
     xlim = c(0, 25), ylim = c(0, 0.3))
lines(density(data_rs_PA$DM_t_PA[which(data_rs_PA$TRT == "NW")]), col = "blue")

# Model initial rosette size at establishment as function of warming treatment
demo_rose <- lmer(DM_t_PA ~ TRT + (1|Group), data = data_rs_PA)

# Perform stepwise selection to minimise AIC
step(demo_rose)

# Minimise AIC by dropping warming term
demo_rose <- lmer(DM_t_PA ~ (1|Group), data = data_rs_PA)
summary(demo_rose)

# Check model assumptions; residuals normally distributed around zero, as expected
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
shapiro.test(resid(demo_rose))
ks.test(resid(demo_rose), pnorm, mean = 0, sd = sd(resid(demo_rose)))
plot(density(resid(demo_rose)), xlim = c(-10, 10), ylim = c(0, 0.25))
plot(demo_rose, xlim = c(8.8, 13.2), ylim = c(-6.5, 6.5))
qqnorm(resid(demo_rose), xlim = c(-3, 3), ylim = c(-7, 7))
qqline(resid(demo_rose))

# Check whether error variance differs significantly between treatment groups; it does not
# Thus, We can model the error term as normal with mean zero and SD agnostic of treatment
plot(density(resid(demo_rose)[which(data_rs_PA$TRT == "W")]), col = "red",
     xlim = c(-10, 10), ylim = c(0, 0.25))
lines(density(resid(demo_rose)[which(data_rs_PA$TRT == "NW")]), col = "blue")
ks.test(resid(demo_rose)[which(data_rs_PA$TRT == "W")],
        resid(demo_rose)[which(data_rs_PA$TRT == "NW")])
var.test(resid(demo_rose)[which(data_rs_PA$TRT == "W")],
         resid(demo_rose)[which(data_rs_PA$TRT == "NW")])

# Store SD of errors to use as stochastic element in demographic simulations
demo_rose_err <- sd(resid(demo_rose))

# Store model coefficients
demo_rose <- fixef(demo_rose)

# Test model fit; seems reasonable
set.seed(979027427)
temp1 <- data_rs_PA$DM_t_PA
temp2 <- rep(demo_rose, length(temp1)) + rnorm(length(temp1), mean = 0, sd = demo_rose_err)
plot(density(temp1), col = "green", xlim = c(0, 25), ylim = c(0, 0.3))
lines(density(temp2), col = "green4")
ks.test(temp1, temp2)

# Remove unused variables
remove(temp1, temp2)





##### Get equations for survival --------------------------------------------------------------------------

# Plot survival versus rosette size; too few deaths for any logistic models to be reliable
# Logistic regression or any sort of MLE would likely lead to small-sample bias
plot(subset(data_rs_PA, TRT == "W")$DM_t_PA, subset(data_rs_PA, TRT == "W")$Survival_PA,
     col = "red", xlim = c(0, 25), ylim = c(0, 1))
points(subset(data_rs_PA, TRT == "NW")$DM_t_PA, subset(data_rs_PA, TRT == "NW")$Survival_PA,
       col = "blue")

# Avoid fit/bias issues by estimating rate as a constant
demo_surv <- mean(data_rs_PA$Survival_PA)





##### Get equations for growth ----------------------------------------------------------------------------

# Plot data; linear model seems like a reasonable choice
temp1 <- data_rs_PA$TRT
temp1[temp1 == "W"] <- "red"
temp1[temp1 == "NW"] <- "blue"
plot(data_rs_PA$DM_t_PA, data_rs_PA$DM_t1_PA,
     xlim = c(0, 25), ylim = c(0, 60), col = temp1)

# Model rosette size at t1 as a function of rosette size at t0; include warming and interaction
demo_grow <- lmer(DM_t1_PA ~ DM_t_PA + TRT + TRT:DM_t_PA + (1|Group),
                  data = data_rs_PA)

# Perform stepwise selection to minimise AIC
step(demo_grow)

# Minimise AIC by dropping interaction term
demo_grow <- lmer(DM_t1_PA ~ DM_t_PA + TRT + (1|Group), data = data_rs_PA)
summary(demo_grow)

# Check model assumptions; residuals normally distributed around zero, as expected
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
shapiro.test(resid(demo_grow))
ks.test(resid(demo_grow), pnorm, mean = 0, sd = sd(resid(demo_grow)))
plot(density(resid(demo_grow)), xlim = c(-15, 15), ylim = c(0, 0.12))
plot(demo_grow, xlim = c(19, 46), ylim = c(-8, 8))
qqnorm(resid(demo_grow), xlim = c(-3, 3), ylim = c(-10, 10))
qqline(resid(demo_grow))

# Check whether error variance differs significantly between treatment groups; it does not
# Thus, We can model the error term as normal with mean zero and SD agnostic of treatment
plot(density(resid(demo_grow)[which(data_rs_PA$TRT == "W")]), col = "red",
     xlim = c(-15, 15), ylim = c(0, 0.12))
lines(density(resid(demo_grow)[which(data_rs_PA$TRT == "NW")]), col = "blue")
ks.test(resid(demo_grow)[which(data_rs_PA$TRT == "W")],
        resid(demo_grow)[which(data_rs_PA$TRT == "NW")])
var.test(resid(demo_grow)[which(data_rs_PA$TRT == "W")],
         resid(demo_grow)[which(data_rs_PA$TRT == "NW")])

# Store SD of errors to use as stochastic element in demographic simulations
demo_grow_err <- sd(resid(demo_grow))

# Store model coefficients
demo_grow_NW <- fixef(demo_grow)[c(1, 2)]
demo_grow_W <- fixef(demo_grow)[c(1, 2)] + c(fixef(demo_grow)[3], 0)

# Test model fit (warmed); seems reasonable
set.seed(386589364)
temp2 <- subset(data_rs_PA, TRT == "W")$DM_t1_PA
temp3 <- demo_grow_W[1] + demo_grow_W[2]*subset(data_rs_PA, TRT == "W")$DM_t_PA +
  rnorm(length(temp2), mean = 0, sd = demo_grow_err)
plot(density(temp2), col = "red", xlim = c(10, 50), ylim = c(0, 0.15))
lines(density(temp3), col = "red4")
ks.test(temp2, temp3)

# Test model fit (unwarmed); seems reasonable
set.seed(386589364)
temp2 <- subset(data_rs_PA, TRT == "NW")$DM_t1_PA
temp3 <- demo_grow_NW[1] + demo_grow_NW[2]*subset(data_rs_PA, TRT == "NW")$DM_t_PA +
  rnorm(length(temp2), mean = 0, sd = demo_grow_err)
plot(density(temp2), col = "blue", xlim = c(10, 50), ylim = c(0, 0.15))
lines(density(temp3), col = "blue4")
ks.test(temp2, temp3)

# Remove unused variables
remove(demo_grow, temp1, temp2, temp3)





##### Get equations for flowering -------------------------------------------------------------------------

# Plot bolting/flowering versus rosette size; too few deaths for any logistic models to be reliable
# Logistic regression or any sort of MLE would likely lead to small-sample bias
plot(subset(data_rs_PA, TRT == "W")$DM_t1_PA, subset(data_rs_PA, TRT == "W")$Flowering_PA,
     col = "red", xlim = c(0, 50), ylim = c(0, 1))
points(subset(data_rs_PA, TRT == "NW")$DM_t1_PA, subset(data_rs_PA, TRT == "NW")$Flowering_PA,
       col = "blue")

# Avoid fit/bias issues by estimating rate as a constant
demo_flow <- mean(data_rs_PA$Flowering_PA)





##### Get equations for number of flower heads ------------------------------------------------------------

# Plot data; point trend is slightly exponential, so log transform response
# Linear model seems like a reasonable choice after transformation
temp1 <- data_rs_PA$TRT
temp1[temp1 == "W"] <- "red"
temp1[temp1 == "NW"] <- "blue"
plot(data_rs_PA$DM_t1_PA, data_rs_PA$Heads_PA,
     xlim = c(0, 50), ylim = c(0, 20), col = temp1)
plot(data_rs_PA$DM_t1_PA, log(data_rs_PA$Heads_PA),
     xlim = c(0, 50), ylim = c(0, 4), col = temp1)

# Model log flower head count as a function of size at t1; include warming and interaction
demo_head <- lmer(log(Heads_PA) ~ DM_t1_PA + TRT + TRT:DM_t1_PA + (1|Group),
                  data = data_rs_PA)

# Perform stepwise selection to minimise AIC
step(demo_head)

# Minimise AIC by dropping interaction term
demo_head <- lmer(log(Heads_PA) ~ DM_t1_PA + TRT + (1|Group), data = data_rs_PA)
summary(demo_head)

# Check model assumptions; residuals normally distributed around zero, as expected
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
shapiro.test(resid(demo_head))
ks.test(resid(demo_head), pnorm, mean = 0, sd = sd(resid(demo_head)))
plot(density(resid(demo_head)), xlim = c(-1.5, 1.5), ylim = c(0, 2))
plot(demo_head, xlim = c(0.8, 3.2), ylim = c(-1.1, 1.1))
qqnorm(resid(demo_head), xlim = c(-3, 3), ylim = c(-1, 1))
qqline(resid(demo_head))

# Check whether error variance differs significantly between treatment groups; it does not
# Thus, We can model the error term as normal with mean zero and SD agnostic of treatment
plot(density(resid(demo_head)[which(data_rs_PA$TRT == "W")]), col = "red",
     xlim = c(-1.5, 1.5), ylim = c(0, 2))
lines(density(resid(demo_head)[which(data_rs_PA$TRT == "NW")]), col = "blue")
ks.test(resid(demo_head)[which(data_rs_PA$TRT == "W")],
        resid(demo_head)[which(data_rs_PA$TRT == "NW")])
var.test(resid(demo_head)[which(data_rs_PA$TRT == "W")],
         resid(demo_head)[which(data_rs_PA$TRT == "NW")])

# Store SD of errors to use as stochastic element in demographic simulations
demo_head_err <- sd(resid(demo_head))

# Store model coefficients
demo_head_NW <- fixef(demo_head)[c(1, 2)]
demo_head_W <- fixef(demo_head)[c(1, 2)] + c(fixef(demo_head)[3], 0)

# Test model fit (warmed); seems reasonable
set.seed(927494733)
temp2 <- subset(data_rs_PA, TRT == "W")$Heads_PA
temp3 <- exp(demo_head_W[1] + demo_head_W[2]*subset(data_rs_PA, TRT == "W")$DM_t1_PA +
               rnorm(length(temp2), mean = 0, sd = demo_head_err))
plot(density(temp2), col = "red", xlim = c(-5, 25), ylim = c(0, 0.3))
lines(density(temp3), col = "red4")
ks.test(temp2, temp3)

# Test model fit (unwarmed); seems reasonable
set.seed(927494733)
temp2 <- subset(data_rs_PA, TRT == "NW")$Heads_PA
temp3 <- exp(demo_head_NW[1] + demo_head_NW[2]*subset(data_rs_PA, TRT == "NW")$DM_t1_PA +
               rnorm(length(temp2), mean = 0, sd = demo_head_err))
plot(density(temp2), col = "blue", xlim = c(-5, 25), ylim = c(0, 0.3))
lines(density(temp3), col = "blue4")
ks.test(temp2, temp3)

# Remove unused variables
remove(demo_head, temp1, temp2, temp3)





##### Get equations for flower head heights ---------------------------------------------------------------

# Plot data; linear model seems like a reasonable choice
temp1 <- data_ht_PA$TRT
temp1[temp1 == "W"] <- "red"
temp1[temp1 == "NW"] <- "blue"
plot(data_ht_PA$DM_t1_PA, data_ht_PA$Height_PA,
     xlim = c(0, 50), ylim = c(0, 145), col = temp1)

# Model flower head height as a function of rosette size; include warming and interaction
# Note the use of t1 size instead of t0; model is thus slightly different than Drees and Shea (2024)
demo_hdht <- lmer(Height_PA ~ DM_t1_PA + TRT + TRT:DM_t1_PA + (1|Group), data = data_ht_PA)

# Perform stepwise selection to minimise AIC
step(demo_hdht)

# Minimise AIC by dropping interaction term
demo_hdht <- lmer(Height_PA ~ DM_t1_PA + TRT + (1|Group), data = data_ht_PA)
summary(demo_hdht)

# Check model assumptions; residuals normally distributed around zero, as expected
# No patterns or heteroskedasticity in residuals also indicates reasonable fit
shapiro.test(resid(demo_hdht))
ks.test(resid(demo_hdht), pnorm, mean = 0, sd = sd(resid(demo_hdht)))
plot(density(resid(demo_hdht)), xlim = c(-40, 40), ylim = c(0, 0.06))
plot(demo_hdht, xlim = c(67, 133), ylim = c(-32, 32))
qqnorm(resid(demo_hdht), xlim = c(-3, 3), ylim = c(-25, 25))
qqline(resid(demo_hdht))

# Estimate errors, but use non-PA distribution since plot averaging mutes small and large heights
# This drastically shrinks the variance of the true distribution of flower head heights
# Underestimating capitulum height variance could mis-estimate dispersal and spread
temp2 <- drop_na(data_ht)
names(temp2)[6] <- c("DM_t1_PA")
temp3 <- predict(demo_hdht, temp2) - temp2$Height

# Check whether error variance differs significantly between treatment groups; it does not
# Thus, We can model the error term as normal with mean zero and SD agnostic of treatment
plot(density(temp3[which(temp2$TRT == "W")]), col = "red",
     xlim = c(-100, 100), ylim = c(0, 0.025))
lines(density(temp3[which(temp2$TRT == "NW")]), col = "blue")
ks.test(temp3[which(temp2$TRT == "W")], temp3[which(temp2$TRT == "NW")])
var.test(temp3[which(temp2$TRT == "W")], temp3[which(temp2$TRT == "NW")])

# Store SD of errors to use as stochastic element in demographic simulations
demo_hdht_err <- sd(temp3)

# Store model coefficients
demo_hdht_NW <- fixef(demo_hdht)[c(1, 2)]
demo_hdht_W <- fixef(demo_hdht)[c(1, 2)] + c(fixef(demo_hdht)[3], 0)

# Test model fit (warmed); seems reasonable
set.seed(284657777)
temp4 <- drop_na(subset(data_ht, TRT == "W"))$Height
temp5 <- demo_hdht_W[1] + demo_hdht_W[2]*drop_na(subset(data_ht, TRT == "W"))$DM_t1 +
  rnorm(length(temp4), mean = 0, sd = demo_hdht_err)
plot(density(temp4), col = "red", xlim = c(0, 200), ylim = c(0, 0.025))
lines(density(temp5), col = "red4")
ks.test(temp4, temp5)

# Test model fit (unwarmed); seems reasonable
set.seed(284657777)
temp6 <- drop_na(subset(data_ht, TRT == "NW"))$Height
temp7 <- demo_hdht_NW[1] + demo_hdht_NW[2]*drop_na(subset(data_ht, TRT == "NW"))$DM_t1 +
  rnorm(length(temp6), mean = 0, sd = demo_hdht_err)
plot(density(temp6), col = "blue", xlim = c(0, 200), ylim = c(0, 0.025))
lines(density(temp7), col = "blue4")
ks.test(temp6, temp7)

# Remove unused variables
remove(demo_hdht, temp1, temp2, temp3, temp4, temp5, temp6, temp7)





##### Get terminal velocity distribution ------------------------------------------------------------------

# Fit lognormal distribution to terminal velocities
# Terminal velocity is drop tube length (1.25 m) divided by drop time
disp_tv_vals <- 1.25/data_tv$drop.time
disp_tv <- fitdistr(disp_tv_vals, "lognormal")$estimate

# Get mean and SD for parameterised lognormal distribution
exp(disp_tv[1] + 0.5*(disp_tv[2]^2))
sqrt(exp(2*disp_tv[1] + (disp_tv[2]^2))*(exp(disp_tv[2]^2) - 1))

# Check that the mean/SD calculations above are accurate; simulated results are close
set.seed(283749842)
mean(rlnorm(1000000, disp_tv[1], disp_tv[2]))
sd(rlnorm(1000000, disp_tv[1], disp_tv[2]))

# Test model fit; seems reasonable
set.seed(283749842)
temp1 <- rlnorm(length(disp_tv_vals), meanlog = disp_tv[1], sdlog = disp_tv[2])
plot(density(disp_tv_vals), col = "green", xlim = c(0, 1.5), ylim = c(0, 4))
lines(density(temp1), col = "green4")
ks.test(disp_tv_vals, temp1)

# Test transformation function, scaling lognormal mean up by 50%; performs as expected
temp2 <- transform.ln(disp_tv[1], disp_tv[2], 1.5, "mean")
exp(disp_tv[1] + 0.5*(disp_tv[2]^2))*1.5
exp(temp2[1] + 0.5*(temp2[2]^2))
sqrt(exp(2*disp_tv[1] + (disp_tv[2]^2))*(exp(disp_tv[2]^2) - 1))
sqrt(exp(2*temp2[1] + (temp2[2]^2))*(exp(temp2[2]^2) - 1))

# Test transformation function, scaling lognormal mean down by 50%; performs as expected
temp3 <- transform.ln(disp_tv[1], disp_tv[2], 0.5, "mean")
exp(disp_tv[1] + 0.5*(disp_tv[2]^2))*0.5
exp(temp3[1] + 0.5*(temp3[2]^2))
sqrt(exp(2*disp_tv[1] + (disp_tv[2]^2))*(exp(disp_tv[2]^2) - 1))
sqrt(exp(2*temp3[1] + (temp3[2]^2))*(exp(temp3[2]^2) - 1))

# Plot mean-scaled distributions against baseline for comparison
set.seed(283749842)
plot(density(rlnorm(1000000, disp_tv[1], disp_tv[2])), xlim = c(0, 2), ylim = c(0, 7))
lines(density(rlnorm(1000000, temp2[1], temp2[2])), col = "green")
lines(density(rlnorm(1000000, temp3[1], temp3[2])), col = "green4")

# Test transformation function, scaling lognormal SD up by 50%; performs as expected
temp4 <- transform.ln(disp_tv[1], disp_tv[2], 1.5, "sd")
exp(disp_tv[1] + 0.5*(disp_tv[2]^2))
exp(temp4[1] + 0.5*(temp4[2]^2))
sqrt(exp(2*disp_tv[1] + (disp_tv[2]^2))*(exp(disp_tv[2]^2) - 1))*1.5
sqrt(exp(2*temp4[1] + (temp4[2]^2))*(exp(temp4[2]^2) - 1))

# Test transformation function, scaling lognormal SD down by 50%; performs as expected
temp5 <- transform.ln(disp_tv[1], disp_tv[2], 0.5, "sd")
exp(disp_tv[1] + 0.5*(disp_tv[2]^2))
exp(temp5[1] + 0.5*(temp5[2]^2))
sqrt(exp(2*disp_tv[1] + (disp_tv[2]^2))*(exp(disp_tv[2]^2) - 1))*0.5
sqrt(exp(2*temp5[1] + (temp5[2]^2))*(exp(temp5[2]^2) - 1))

# Plot SD-scaled distributions against baseline for comparison
set.seed(283749842)
plot(density(rlnorm(1000000, disp_tv[1], disp_tv[2])), xlim = c(0, 2), ylim = c(0, 7))
lines(density(rlnorm(1000000, temp4[1], temp4[2])), col = "green")
lines(density(rlnorm(1000000, temp5[1], temp5[2])), col = "green4")

# Remove unused variables
remove(temp1, temp2, temp3, temp4, temp5)





##### Get wind speed distribution -------------------------------------------------------------------------

# Fit Weibull distribution to wind speeds
# Assume no seed release occurs for wind speeds of zero, so remove zero values
disp_ws_vals <- c(data_ws1$Wind1, data_ws2$Wind1)
disp_ws_vals <- disp_ws_vals[disp_ws_vals > 0]
disp_ws <- fitdistr(disp_ws_vals, "weibull")$estimate

# Get mean and SD for parameterised Weibull distribution
disp_ws[2]*gamma(1 + 1/disp_ws[1])
sqrt((disp_ws[2]^2)*(gamma(1 + 2/disp_ws[1]) - gamma(1 + 1/disp_ws[1])^2))

# Check that the mean/SD calculations above are accurate; simulated results are close
set.seed(749402744)
mean(rweibull(1000000, disp_ws[1], disp_ws[2]))
sd(rweibull(1000000, disp_ws[1], disp_ws[2]))

# Test model fit; seems reasonable
# Wind data is a bit messy and sample size is very high, so K-S is not significant
# This is fine, though, as many studies have used Weibull for wind speed distribution
set.seed(749402744)
temp1 <- rweibull(length(disp_ws_vals), disp_ws[1], disp_ws[2])
plot(density(disp_ws_vals), col = "green", xlim = c(-1, 15), ylim = c(0, 0.4))
lines(density(temp1), col = "green4")
ks.test(disp_ws_vals, temp1)

# Test transformation function, scaling Weibull mean up by 50%; performs as expected
temp2 <- transform.wb(disp_ws[1], disp_ws[2], 1.5, "mean")
disp_ws[2]*gamma(1 + 1/disp_ws[1])*1.5
temp2[2]*gamma(1 + 1/temp2[1])
sqrt((disp_ws[2]^2)*(gamma(1 + 2/disp_ws[1]) - gamma(1 + 1/disp_ws[1])^2))
sqrt((temp2[2]^2)*(gamma(1 + 2/temp2[1]) - gamma(1 + 1/temp2[1])^2))

# Test transformation function, scaling Weibull mean down by 50%; performs as expected
temp3 <- transform.wb(disp_ws[1], disp_ws[2], 0.5, "mean")
disp_ws[2]*gamma(1 + 1/disp_ws[1])*0.5
temp3[2]*gamma(1 + 1/temp3[1])
sqrt((disp_ws[2]^2)*(gamma(1 + 2/disp_ws[1]) - gamma(1 + 1/disp_ws[1])^2))
sqrt((temp3[2]^2)*(gamma(1 + 2/temp3[1]) - gamma(1 + 1/temp3[1])^2))

# Plot mean-scaled distributions against baseline for comparison
set.seed(749402744)
plot(density(rweibull(1000000, disp_ws[1], disp_ws[2])), xlim = c(0, 15), ylim = c(0, 1))
lines(density(rweibull(1000000, temp2[1], temp2[2])), col = "green")
lines(density(rweibull(1000000, temp3[1], temp3[2])), col = "green4")

# Test transformation function, scaling Weibull SD up by 50%; performs as expected
temp4 <- transform.wb(disp_ws[1], disp_ws[2], 1.5, "sd")
disp_ws[2]*gamma(1 + 1/disp_ws[1])
temp4[2]*gamma(1 + 1/temp4[1])
sqrt((disp_ws[2]^2)*(gamma(1 + 2/disp_ws[1]) - gamma(1 + 1/disp_ws[1])^2))*1.5
sqrt((temp4[2]^2)*(gamma(1 + 2/temp4[1]) - gamma(1 + 1/temp4[1])^2))

# Test transformation function, scaling Weibull SD down by 50%; performs as expected
temp5 <- transform.wb(disp_ws[1], disp_ws[2], 0.5, "sd")
disp_ws[2]*gamma(1 + 1/disp_ws[1])
temp5[2]*gamma(1 + 1/temp5[1])
sqrt((disp_ws[2]^2)*(gamma(1 + 2/disp_ws[1]) - gamma(1 + 1/disp_ws[1])^2))*0.5
sqrt((temp5[2]^2)*(gamma(1 + 2/temp5[1]) - gamma(1 + 1/temp5[1])^2))

# Plot SD-scaled distributions against baseline for comparison
set.seed(749402744)
plot(density(rweibull(1000000, disp_ws[1], disp_ws[2])), xlim = c(0, 15), ylim = c(0, 1))
lines(density(rweibull(1000000, temp4[1], temp4[2])), col = "green")
lines(density(rweibull(1000000, temp5[1], temp5[2])), col = "green4")

# Remove unused variables
remove(temp1, temp2, temp3, temp4, temp5)

