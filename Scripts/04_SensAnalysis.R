##### Run simulations for base cases ----------------------------------------------------------------------

# Load simulation settings before running any wavespeed code
sets_plot <- FALSE     # (De)Activate wave plotting
sets_nest <- TRUE      # (De)Activate ant nests
sets_trim <- TRUE      # (De)Activate core area trimming
sets_years <- 1000     # Number of years to simulate
sets_trimD <- 1000     # Distance (m) behind wavefront to trim

# Set demography/dispersal parameters not covered in ADSP or WDSP
misc_tDens <- 10       # Max thistle density per metre
misc_nDens <- 0.2      # Ant nest density (nests/m)
misc_nRange <- 25      # Max detection range (m) from ant nests

# Set up ADSP parameter transformation for runs unwarmed, no ant dispersal or predation
sens_w_adsp_num1 <- c(20, 21)
sens_w_adsp_val1 <- c(0, 10)

# Set up ADSP parameter transformation for runs warmed, no ant dispersal or predation
sens_w_adsp_num2 <- c(4, 8, 11, 16, 18, 20, 21)
sens_w_adsp_val2 <- c(demo_grow_W[1]/demo_grow_NW[1], demo_head_W[1]/demo_head_NW[1],
                      demo_hdht_W[1]/demo_hdht_NW[1], 0.302/0.233, 0.302/0.233, 0, 10)

# Set up ADSP parameter transformation for runs warmed, ant dispersal and predation
sens_w_adsp_num3 <- c(4, 8, 11, 16, 18, 20)
sens_w_adsp_val3 <- c(demo_grow_W[1]/demo_grow_NW[1], demo_head_W[1]/demo_head_NW[1],
                      demo_hdht_W[1]/demo_hdht_NW[1], 0.302/0.233, 0.302/0.233, 0.984/0.948)

# Set up WDSP parameter transformation runs warmed, ant dispersal and predation irrelevant
sens_w_wdsp_num1 <- c(1)
sens_w_wdsp_val1 <- c(0.129/0.056)

# Note: for each of the code chunks below, run 03_Wavepseeds after setting parameters
# Then output wavefront data to CSV for later analysis

# [1] Base case: unwarmed, no ant dispersal or predation
set.seed(23568)
sets_nest <- FALSE
aVec <- adsp.param(aNum = sens_w_adsp_num1, aVal = sens_w_adsp_val1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_base_1.csv")

# [2] Base case: warmed, no ant dispersal or predation
set.seed(72363)
sets_nest <- FALSE
aVec <- adsp.param(aNum = sens_w_adsp_num2, aVal = sens_w_adsp_val2)
wVec <- wdsp.param(wNum = sens_w_wdsp_num1, wVal = sens_w_wdsp_val1)
write.csv(wv_front, "Data/wv_base_2.csv")

# [3] Base case: unwarmed, ant dispersal and predation
set.seed(92362)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_base_3.csv")

# [4] Base case: warmed, ant dispersal and predation
set.seed(83574)
aVec <- adsp.param(aNum = sens_w_adsp_num3, aVal = sens_w_adsp_val3)
wVec <- wdsp.param(wNum = sens_w_wdsp_num1, wVal = sens_w_wdsp_val1)
write.csv(wv_front, "Data/wv_base_4.csv")





##### Run simulations for various cases across dispersal parameter space ----------------------------------

# Load simulation settings before running any wavespeed code
sets_plot <- FALSE     # (De)Activate wave plotting
sets_nest <- TRUE      # (De)Activate ant nests
sets_trim <- TRUE      # (De)Activate core area trimming
sets_years <- 100      # Number of years to simulate
sets_trimD <- 1000     # Distance (m) behind wavefront to trim

# Set demography/dispersal parameters not covered in ADSP or WDSP
misc_tDens <- 10       # Max thistle density per metre
misc_nDens <- 0.2      # Ant nest density (nests/m)
misc_nRange <- 25      # Max detection range (m) from ant nests

# Set span of multipliers for parameters not strictly in space [0, 1]
# Set span of probabilities/rates for parameters strictly in space [0, 1]
multPr <- seq(0.0, 1.0, 0.1)
multRe <- seq(0.5, 1.5, 0.1)

# [1] Modify head height intercept
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- lapply(X = multRe, FUN = adsp.param, aNum = 11)
write.csv(wv_span, "Data/wv_span_1.csv")

# [2] Modify head height size-slope
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- lapply(X = multRe, FUN = adsp.param, aNum = 12)
write.csv(wv_span, "Data/wv_span_2.csv")

# [3] Modify head height SD
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- lapply(X = multRe, FUN = adsp.param, aNum = 13)
write.csv(wv_span, "Data/wv_span_3.csv")

# [4] Modify prob. of seed removal by ants
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- lapply(X = multPr, FUN = adsp.param, aNum = 20)
write.csv(wv_span, "Data/wv_span_4.csv")

# [5] Modify prob. of surviving predation if removed by ants
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- lapply(X = multPr, FUN = adsp.param, aNum = 21)
write.csv(wv_span, "Data/wv_span_5.csv")

# [6] Modify prob. of seed release from capitulum
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multPr, FUN = wdsp.param, wNum = 1)
write.csv(wv_span, "Data/wv_span_6.csv")

# [7] Modify surrounding vegetation height
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multRe, FUN = wdsp.param, wNum = 2)
write.csv(wv_span, "Data/wv_span_7.csv")

# [8] Modify mean seed terminal velocity
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multRe, FUN = wdsp.param, wNum = 3)
write.csv(wv_span, "Data/wv_span_8.csv")

# [9] Modify SD seed terminal velocity
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multRe, FUN = wdsp.param, wNum = 4)
write.csv(wv_span, "Data/wv_span_9.csv")

# [10] Modify mean wind speed
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multRe, FUN = wdsp.param, wNum = 5)
write.csv(wv_span, "Data/wv_span_10.csv")

# [11] Modify SD wind speed
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multRe, FUN = wdsp.param, wNum = 6)
write.csv(wv_span, "Data/wv_span_11.csv")

# [12] Modify max thistle density
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- 10*multRe
write.csv(wv_span, "Data/wv_span_12.csv")

# [13] Modify ant nest density 
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
WVecList <- 0.2*multRe
write.csv(wv_span, "Data/wv_span_13.csv")

# [14] Modify max detection range from ant nests
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
WVecList <- 25*multRe
write.csv(wv_span, "Data/wv_span_14.csv")

