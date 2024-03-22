##### Base runs -------------------------------------------------------------------------------------------

# Ensure parameters below are set before running any of options 1-4
nestOn <- TRUE    # Should ant dispersal be included?
trim <- TRUE      # Should core area of wave be trimmed?
plotOn <- FALSE   # Should the wave be plotted?
nYear <- 1000     # Number of years to simulate
trimAmt <- 1000   # Distance (m) behind wavefront to trim
tDens <- 10       # Max thistle density per metre
range <- 25       # Max detection range (m) from ant nests
nDens <- 0.2      # Ant nest density (nests/m)

# If running a warmed scenario, run the code below to update demographic/dispersal parameters
# Asterisks indicate parameters that are different under warming conditions
wdsp.param <- function(wNum, wVal){
  wParam <- c()
  wParam[1] <- 0.15                 # Vegetation height (m)
  wParam[2] <- disp_ws[1]           # Shape wind speed, Weibull dist.
  wParam[3] <- disp_ws[2]           # Scale wind speed, Weibull dist.
  wParam[4] <- disp_tv[1]           # Log mean terminal velocity, lognormal dist.
  wParam[5] <- disp_tv[2]           # Log SD terminal velocity, lognormal dist.
  wParam[6] <- 0.129                # Probability of seed release from capitulum***
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
  return(wParam)}
adsp.param <- function(aNum, aVal){
  aParam <- c()
  aParam[1] <- 476                  # Seeds per flower head
  aParam[2] <- 0.850                # Prob. of surviving pre-dispersal seed predation (florivory)
  aParam[3] <- 0.984                # Prob. of seed removal by ants***
  aParam[4] <- 0.100                # Prob. of surviving post-dispersal seed predation by ants
  aParam[5] <- 0.302                # Prob. of establishment from seed***
  aParam[6] <- 0.233                # Prob. of seed entering seed bank
  aParam[7] <- 0.302                # Prob. of seed establishing from seed bank***
  aParam[8] <- 0.260                # Prob. of seed survival in seed bank
  aParam[9] <- demo_rose            # Rosette size t0 intercept
  aParam[10] <- demo_rose_err       # Rosette size t0 model SD
  aParam[11] <- demo_surv           # Prob. survival
  aParam[12] <- demo_grow_W[1]      # Rosette size t1 intercept***
  aParam[13] <- demo_grow_W[2]      # Rosette size t1 size-slope
  aParam[14] <- demo_grow_err       # Mean growth model SD
  aParam[15] <- demo_flow           # Prob. flowering
  aParam[16] <- demo_head_W[1]      # Num. heads intercept***
  aParam[17] <- demo_head_W[2]      # Num. heads size-slope
  aParam[18] <- demo_head_err       # Num. heads model SD
  aParam[19] <- demo_hdht_W[1]      # Head height intercept***
  aParam[20] <- demo_hdht_W[2]      # Head height diameter-slope
  aParam[21] <- demo_hdht_err       # Head height model SD
  aVec <- rep(1, length(aParam))
  aVec[aNum] <- aVal
  aParam <- aParam*aVec
  return(aParam)}

# [1] Seed predation enabled, ant dispersal enabled [unwarmed]
set.seed(28360)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)

# [2] Seed predation enabled, ant dispersal enabled [warmed]
set.seed(58334)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)

# [3] Seed predation disabled, ant dispersal disabled [unwarmed]
set.seed(80011)
aVec <- adsp.param(aNum = 4, aVal = 10)
wVec <- wdsp.param(wNum = 1, wVal = 1)
nestOn <- FALSE

# [4] Seed predation disabled, ant dispersal disabled [unwarmed]
set.seed(42842)
aVec <- adsp.param(aNum = 4, aVal = 10)
wVec <- wdsp.param(wNum = 1, wVal = 1)
nestOn <- FALSE





##### Sensitivity analyses --------------------------------------------------------------------------------

# Ensure parameters below are set before running any of options 1-30
nestOn <- TRUE    # Should ant dispersal be included?
trim <- TRUE      # Should core area of wave be trimmed?
plotOn <- FALSE   # Should the wave be plotted?
nYear <- 1000     # Number of years to simulate
trimAmt <- 1000   # Distance (m) behind wavefront to trim
tDens <- 10       # Max thistle density per metre
range <- 25       # Max detection range (m) from ant nests
nDens <- 0.2      # Ant nest density (nests/m)

# [1] Seeds per flower head
set.seed(63925)
aVec <- adsp.param(aNum = 1, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_1.csv")
write.csv(wv_dens, "Data/wv_dens_1.csv")

# [2] Prob. of surviving pre-dispersal seed predation (florivory)
set.seed(12263)
aVec <- adsp.param(aNum = 2, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_2.csv")
write.csv(wv_dens, "Data/wv_dens_2.csv")

# [3] Prob. of seed removal by ants
set.seed(98374)
aVec <- adsp.param(aNum = 3, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_3.csv")
write.csv(wv_dens, "Data/wv_dens_3.csv")

# [4] Prob. of surviving post-dispersal seed predation by ants
set.seed(29342)
aVec <- adsp.param(aNum = 4, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_4.csv")
write.csv(wv_dens, "Data/wv_dens_4.csv")

# [5] Prob. of establishment from seed
set.seed(36543)
aVec <- adsp.param(aNum = 5, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_5.csv")
write.csv(wv_dens, "Data/wv_dens_5.csv")

# [6] Prob. of seed entering seed bank
set.seed(14563)
aVec <- adsp.param(aNum = 6, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_6.csv")
write.csv(wv_dens, "Data/wv_dens_6.csv")

# [7] Prob. of seed establishing from seed bank
set.seed(77253)
aVec <- adsp.param(aNum = 7, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_7.csv")
write.csv(wv_dens, "Data/wv_dens_7.csv")

# [8] Prob. of seed survival in seed bank
set.seed(25657)
aVec <- adsp.param(aNum = 8, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_8.csv")
write.csv(wv_dens, "Data/wv_dens_8.csv")

# [9] Rosette size t0 intercept
set.seed(87324)
aVec <- adsp.param(aNum = 9, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_9.csv")
write.csv(wv_dens, "Data/wv_dens_9.csv")

# [10] Rosette size t0 model SD
set.seed(55695)
aVec <- adsp.param(aNum = 10, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_10.csv")
write.csv(wv_dens, "Data/wv_dens_10.csv")

# [11] Prob. survival
set.seed(63453)
aVec <- adsp.param(aNum = 11, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_11.csv")
write.csv(wv_dens, "Data/wv_dens_11.csv")

# [12] Rosette size t1 intercept
set.seed(45675)
aVec <- adsp.param(aNum = 12, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_12.csv")
write.csv(wv_dens, "Data/wv_dens_12.csv")

# [13] Rosette size t1 size-slope
set.seed(34675)
aVec <- adsp.param(aNum = 13, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_13.csv")
write.csv(wv_dens, "Data/wv_dens_13.csv")

# [14] Mean growth model SD
set.seed(90344)
aVec <- adsp.param(aNum = 14, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_14.csv")
write.csv(wv_dens, "Data/wv_dens_14.csv")

# [15] Prob. flowering
set.seed(47333)
aVec <- adsp.param(aNum = 15, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_15.csv")
write.csv(wv_dens, "Data/wv_dens_15.csv")

# [16] Num. heads intercept
set.seed(69453)
aVec <- adsp.param(aNum = 16, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_16.csv")
write.csv(wv_dens, "Data/wv_dens_16.csv")

# [17] Num. heads size-slope
set.seed(55484)
aVec <- adsp.param(aNum = 17, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_17.csv")
write.csv(wv_dens, "Data/wv_dens_17.csv")

# [18] Num. heads model SD
set.seed(73454)
aVec <- adsp.param(aNum = 18, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_18.csv")
write.csv(wv_dens, "Data/wv_dens_18.csv")

# [19] Head height intercept
set.seed(25955)
aVec <- adsp.param(aNum = 19, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_19.csv")
write.csv(wv_dens, "Data/wv_dens_19.csv")

# [20] Head height diameter-slope
set.seed(14686)
aVec <- adsp.param(aNum = 20, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_20.csv")
write.csv(wv_dens, "Data/wv_dens_20.csv")

# [21] Head height model SD
set.seed(82450)
aVec <- adsp.param(aNum = 21, aVal = 0.9)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_front_21.csv")
write.csv(wv_dens, "Data/wv_dens_21.csv")

# [22] Vegetation height (m)
set.seed(14578)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 0.9)
write.csv(wv_front, "Data/wv_front_22.csv")
write.csv(wv_dens, "Data/wv_dens_22.csv")

# [23] Mean wind speed, Weibull dist.
set.seed(51487)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 2, wVal = 0.9)
write.csv(wv_front, "Data/wv_front_23.csv")
write.csv(wv_dens, "Data/wv_dens_23.csv")

# [24] SD wind speed, Weibull dist.
set.seed(93106)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 3, wVal = 0.9)
write.csv(wv_front, "Data/wv_front_24.csv")
write.csv(wv_dens, "Data/wv_dens_24.csv")

# [25] Mean terminal velocity, lognormal dist.
set.seed(31903)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 4, wVal = 0.9)
write.csv(wv_front, "Data/wv_front_25.csv")
write.csv(wv_dens, "Data/wv_dens_25.csv")

# [26] SD terminal velocity, lognormal dist.
set.seed(72074)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 5, wVal = 0.9)
write.csv(wv_front, "Data/wv_front_26.csv")
write.csv(wv_dens, "Data/wv_dens_26.csv")

# [27] Probability of seed release from capitulum
set.seed(49921)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 6, wVal = 0.9)
write.csv(wv_front, "Data/wv_front_27.csv")
write.csv(wv_dens, "Data/wv_dens_27.csv")

# [28] Max thistle density per metre
set.seed(52740)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
tDens <- 10 
write.csv(wv_front, "Data/wv_front_28.csv")
write.csv(wv_dens, "Data/wv_dens_28.csv")

# [29] Max detection range from ant nests
set.seed(45855)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
range <- 25
write.csv(wv_front, "Data/wv_front_29.csv")
write.csv(wv_dens, "Data/wv_dens_29.csv")

# [30] Ant nest density
set.seed(12568)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
nDens <- 0.2
write.csv(wv_front, "Data/wv_front_30.csv")
write.csv(wv_dens, "Data/wv_dens_30.csv")

