##### Base runs -------------------------------------------------------------------------------------------

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

# Set up ADSP parameter transformation for runs [unwarmed, no ant dispersal/predation]
sens_w_adsp_num1 <- c(20, 21)
sens_w_adsp_val1 <- c(0, 10)

# Set up ADSP parameter transformation for runs [warmed, no ant dispersal/predation]
sens_w_adsp_num2 <- c(4, 8, 11, 16, 18, 20, 21)
sens_w_adsp_val2 <- c(demo_grow_W[1]/demo_grow_NW[1], demo_head_W[1]/demo_head_NW[1],
                      demo_hdht_W[1]/demo_hdht_NW[1], 0.302/0.233, 0.302/0.233, 0, 10)

# Set up ADSP parameter transformation for runs [warmed, ant dispersal/predation]
sens_w_adsp_num3 <- c(4, 8, 11, 16, 18, 20)
sens_w_adsp_val3 <- c(demo_grow_W[1]/demo_grow_NW[1], demo_head_W[1]/demo_head_NW[1],
                      demo_hdht_W[1]/demo_hdht_NW[1], 0.302/0.233, 0.302/0.233, 0.984/0.948)

# Set up WDSP parameter transformation runs [warmed]
sens_w_wdsp_num1 <- c(1)
sens_w_wdsp_val1 <- c(0.129/0.056)

# [1] Base run with no ant dispersal or predation [unwarmed]
set.seed()
sets_nest <- FALSE
aVec <- adsp.param(aNum = sens_w_adsp_num1, aVal = sens_w_adsp_val1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_base_1.csv")

# [2] Base run with no ant dispersal or predation [warmed]
set.seed()
sets_nest <- FALSE
aVec <- adsp.param(aNum = sens_w_adsp_num2, aVal = sens_w_adsp_val2)
wVec <- wdsp.param(wNum = sens_w_wdsp_num1, wVal = sens_w_wdsp_val1)
write.csv(wv_front, "Data/wv_base_2.csv")

# [3] Base run with ant dispersal and predation [unwarmed]
set.seed()
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
write.csv(wv_front, "Data/wv_base_3.csv")

# [4] Base run with ant dispersal and predation [warmed]
set.seed()
aVec <- adsp.param(aNum = sens_w_adsp_num3, aVal = sens_w_adsp_val3)
wVec <- wdsp.param(wNum = sens_w_wdsp_num1, wVal = sens_w_wdsp_val1)
write.csv(wv_front, "Data/wv_base_4.csv")





##### Parameter exploration runs --------------------------------------------------------------------------

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

# Set multipliers and baseline parameters
multPr <- seq(0.0, 1.0, 0.1)
multRe <- seq(0.5, 1.5, 0.1)

# [1] Head height intercept
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- lapply(X = multRe, FUN = adsp.param, aNum = 11)
write.csv(wv_span, "Data/wv_span_2.csv")

# [2] Head height size-slope
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- lapply(X = multRe, FUN = adsp.param, aNum = 12)
write.csv(wv_span, "Data/wv_span_3.csv")

# [3] Head height size-slope
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- lapply(X = multRe, FUN = adsp.param, aNum = 13)
write.csv(wv_span, "Data/wv_span_4.csv")

# [4] Prob. of seed removal by ants
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- lapply(X = multPr, FUN = adsp.param, aNum = 20)
write.csv(wv_span, "Data/wv_span_1.csv")

# [5] Prob. of surviving predation if removed by ants
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- lapply(X = multPr, FUN = adsp.param, aNum = 21)
write.csv(wv_span, "Data/wv_span_2.csv")

# [6] Prob. of seed release from capitulum
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multPr, FUN = wdsp.param, wNum = 1)
write.csv(wv_span, "Data/wv_span_2.csv")

# [7] Surrounding vegetation height
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multRe, FUN = wdsp.param, wNum = 2)
write.csv(wv_span, "Data/wv_span_2.csv")

# [8] Mean seed terminal velocity
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multRe, FUN = wdsp.param, wNum = 3)
write.csv(wv_span, "Data/wv_span_2.csv")

# [9] SD seed terminal velocity
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multRe, FUN = wdsp.param, wNum = 4)
write.csv(wv_span, "Data/wv_span_2.csv")

# [10] Mean wind speed
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multRe, FUN = wdsp.param, wNum = 5)
write.csv(wv_span, "Data/wv_span_2.csv")

# [11] SD wind speed
aVec <- adsp.param(aNum = 1, aVal = 1)
simSeeds <- c()
wVecList <- lapply(X = multRe, FUN = wdsp.param, wNum = 6)
write.csv(wv_span, "Data/wv_span_2.csv")

# [12] Max thistle density
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
aVecList <- 10*multRe

# [13] Ant nest density 
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
WVecList <- 0.2*multRe

# [14] Max detection range from ant nests
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
simSeeds <- c()
WVecList <- 25*multRe
