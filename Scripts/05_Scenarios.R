##### Run simulations for base cases ----------------------------------------------------------------------

# Load simulation settings before running any wavespeed code
sets_plot <- FALSE     # (De)Activate wave plotting
sets_trim <- TRUE      # (De)Activate core area trimming
sets_years <- 1000     # Number of years to simulate
sets_trimD <- 1000     # Distance (m) behind wavefront to trim
sets_dCnts <- TRUE     # Track wind/ant dispersal contributions for each seed

# Set demography/dispersal parameters not covered in ADSP or WDSP
misc_tDens <- 10       # Max thistle density per metre
misc_nDens <- 0.2      # Ant nest density (nests/m)
misc_nRange <- 25      # Max detection range (m) from ant nests

# Set up ADSP parameter transformation: warmed w/ no ant predation
sens_w_adsp_num1 <- c(4, 8, 11, 16, 18, 20, 21)
sens_w_adsp_val1 <- c(demo_grow_W[1]/demo_grow_NW[1], demo_head_W[1]/demo_head_NW[1],
                      demo_hdht_W[1]/demo_hdht_NW[1], 0.302/0.233, 0.302/0.233, 0.984/0.948, 0)

# Set up ADSP parameter transformation: warmed w/ ant predation
sens_w_adsp_num2 <- c(4, 8, 11, 16, 18, 20)
sens_w_adsp_val2 <- c(demo_grow_W[1]/demo_grow_NW[1], demo_head_W[1]/demo_head_NW[1],
                      demo_hdht_W[1]/demo_hdht_NW[1], 0.302/0.233, 0.302/0.233, 0.984/0.948)

# For each code chunk below, first run the lines from seed set to wVec definition (inclusive)
# Then run the entirety of 03_Wavespeeds and calculate median wavespeed
# Optional: save wavespeeds to CSV so there's no need to re-run code

# [1] Base case: unwarmed w/ no ant dispersal or predation
# Note: reset misc_nDens and misc_nRange before running any other code blocks
set.seed(68932)
misc_nDens <- 0  
misc_nRange <- 0
aVec <- adsp.param(aNum = 21, aVal = 0)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_base_1.csv")
write.csv(wv_dCnts, "Data/wv_base_1d.csv")
write.csv(wv_disp_w, "Data/wv_base_1w.csv")
write.csv(wv_disp_a, "Data/wv_base_1a.csv")

# [2] Base case: warmed w/ no ant dispersal or predation
# Note: reset misc_nDens and misc_nRange before running any other code blocks
set.seed(72363)
misc_nDens <- 0
misc_nRange <- 0
aVec <- adsp.param(aNum = sens_w_adsp_num1, aVal = sens_w_adsp_val1)
wVec <- wdsp.param(wNum = 1, wVal = 0.129/0.056)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_base_2.csv")
write.csv(wv_dCnts, "Data/wv_base_2d.csv")
write.csv(wv_disp_w, "Data/wv_base_2w.csv")
write.csv(wv_disp_a, "Data/wv_base_2a.csv")

# [3] Base case: unwarmed w/ ant dispersal and predation
set.seed(10845)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_base_3.csv")
write.csv(wv_dCnts, "Data/wv_base_3d.csv")
write.csv(wv_disp_w, "Data/wv_base_3w.csv")
write.csv(wv_disp_a, "Data/wv_base_3a.csv")

# [4] Base case: warmed w/ ant dispersal and predation
set.seed(92376)
aVec <- adsp.param(aNum = sens_w_adsp_num2, aVal = sens_w_adsp_val2)
wVec <- wdsp.param(wNum = 1, wVal = 0.129/0.056)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_base_4.csv")
write.csv(wv_dCnts, "Data/wv_base_4d.csv")
write.csv(wv_disp_w, "Data/wv_base_4w.csv")
write.csv(wv_disp_a, "Data/wv_base_4a.csv")

# [5] Base case: unwarmed w/ no ant dispersal, but predation
# Note: reset misc_nDens and misc_nRange before running any other code blocks
set.seed(68253)
misc_nDens <- 0
misc_nRange <- 0 
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_base_5.csv")
write.csv(wv_dCnts, "Data/wv_base_5d.csv")
write.csv(wv_disp_w, "Data/wv_base_5w.csv")
write.csv(wv_disp_a, "Data/wv_base_5a.csv")

# [6] Base case: warmed w/ no ant dispersal, but predation
# Note: reset misc_nDens and misc_nRange before running any other code blocks
set.seed(32674)
misc_nDens <- 0
misc_nRange <- 0 
aVec <- adsp.param(aNum = sens_w_adsp_num2, aVal = sens_w_adsp_val2)
wVec <- wdsp.param(wNum = 1, wVal = 0.129/0.056)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_base_6.csv")
write.csv(wv_dCnts, "Data/wv_base_6d.csv")
write.csv(wv_disp_w, "Data/wv_base_6w.csv")
write.csv(wv_disp_a, "Data/wv_base_6a.csv")

# [7] Base case: unwarmed w/ ant dispersal, but no predation
set.seed(93265)
aVec <- adsp.param(aNum = 21, aVal = 0)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_base_7.csv")
write.csv(wv_dCnts, "Data/wv_base_7d.csv")
write.csv(wv_disp_w, "Data/wv_base_7w.csv")
write.csv(wv_disp_a, "Data/wv_base_7a.csv")

# [8] Base case: warmed w/ ant dispersal, but no predation
set.seed(40911)
aVec <- adsp.param(aNum = sens_w_adsp_num1, aVal = sens_w_adsp_val1)
wVec <- wdsp.param(wNum = 1, wVal = 0.129/0.056)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_base_8.csv")
write.csv(wv_dCnts, "Data/wv_base_8d.csv")
write.csv(wv_disp_w, "Data/wv_base_8w.csv")
write.csv(wv_disp_a, "Data/wv_base_8a.csv")

# [9] Base case: unwarmed w/ ant dispersal, but no predation [1 nest/m]
# Shows that density effects are what's driving decreased speeds when ant dispersal is active
# Note: reset misc_nDens before running any other code blocks
set.seed(38563)
misc_nDens <- 1
aVec <- adsp.param(aNum = 21, aVal = 0)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_base_9.csv")

# [10] Base case: unwarmed w/ ant dispersal and predation [1 nest/m]
# Shows that density effects are what's driving decreased speeds when ant dispersal is active
# Note: reset misc_nDens before running any other code blocks
set.seed(91865)
misc_nDens <- 1
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_base_10.csv")

# Load in manually-combined version of CSVs above for future use
# Includes only median wavespeeds, not dispersal stats
wv_base <- read.csv("Data/wv_base.csv")





##### Run simulations for perturbation analysis -----------------------------------------------------------

# Load simulation settings before running any wavespeed code
sets_plot <- FALSE     # (De)Activate wave plotting
sets_trim <- TRUE      # (De)Activate core area trimming
sets_years <- 100      # Number of years to simulate
sets_trimD <- 1000     # Distance (m) behind wavefront to trim
sets_dCnts <- FALSE    # Track wind/ant dispersal contributions for each seed

# Set demography/dispersal parameters not covered in ADSP or WDSP
misc_tDens <- 10       # Max thistle density per metre
misc_nDens <- 0.2      # Ant nest density (nests/m)
misc_nRange <- 25      # Max detection range (m) from ant nests

# Set span of probabilities/rates for parameters strictly in space [0, 1]
# Set span of multipliers for parameters not strictly in space [0, 1]
pert_prob <- seq(0.0, 1.0, 0.1)
pert_real <- seq(0.5, 1.5, 0.1)

# For each code chunk below, first run the lines from pert_seed to wv_pert definitions (inclusive)
# Then perturb parameter in question; see notes at bottom of section for details
# Next, set the first seed using set.seed(pert_seed[1])
# Then run the entirety of 03_Wavespeeds for the first rep
# Next, calculate median wavespeed and append to wv_pert
# Repeat parameter perturbation, seed set, and script run for reps 2-11
# Optional: save wavespeeds to CSV so there's no need to re-run code

# [1] Modify head height intercept
pert_seed <- c(36173, 48342, 03421, 11734, 93615, 12362, 77715, 54161, 89235, 91026, 39973)
wVec <- wdsp.param(wNum = 1, wVal = 1)
aVecList <- lapply(X = pert_real, FUN = adsp.param, aNum = 11)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_1.csv")

# [2] Modify head height size-slope
pert_seed <- c(01974, 82353, 91865, 18765, 19651, 29374, 71624, 65192, 51256, 83763, 28262)
wVec <- wdsp.param(wNum = 1, wVal = 1)
aVecList <- lapply(X = pert_real, FUN = adsp.param, aNum = 12)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_2.csv")

# [3] Modify head height SD
pert_seed <- c(23522, 83563, 10945, 12723, 28936, 06256, 29753, 00583, 49525, 20017, 10764)
wVec <- wdsp.param(wNum = 1, wVal = 1)
aVecList <- lapply(X = pert_real, FUN = adsp.param, aNum = 13)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_3.csv")

# [4] Modify prob. of seed removal by ants
pert_seed <- c(89345, 19015, 79263, 72957, 81253, 61953, 13533, 91566, 49176, 19673, 20194)
wVec <- wdsp.param(wNum = 1, wVal = 1)
aVecList <- lapply(X = pert_prob/0.948, FUN = adsp.param, aNum = 20)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_4.csv")

# [5] Modify prob. of predation if removed by ants
pert_seed <- c(68385, 91064, 07354, 62365, 20856, 89264, 91907, 60175, 10456, 29013, 29222)
wVec <- wdsp.param(wNum = 1, wVal = 1)
aVecList <- lapply(X = pert_prob/0.900, FUN = adsp.param, aNum = 21)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_5.csv")

# [6] Modify prob. of seed release from capitulum
pert_seed <- c(93664, 08523, 18526, 38566, 59725, 82965, 99365, 77724, 28962, 20911, 92651)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVecList <- lapply(X = pert_prob/0.056, FUN = wdsp.param, wNum = 1)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_6.csv")

# [7] Modify surrounding vegetation height
pert_seed <- c(60923, 49622, 58915, 09275, 40825, 90141, 80000, 76252, 91964, 69862, 27151)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVecList <- lapply(X = pert_real, FUN = wdsp.param, wNum = 2)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_7.csv")

# [8] Modify mean seed terminal velocity
pert_seed <- c(88325, 29376, 93646, 83577, 72634, 16511, 40174, 50285, 39754, 33753, 58333)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVecList <- lapply(X = pert_real, FUN = wdsp.param, wNum = 3)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_8.csv")

# [9] Modify SD seed terminal velocity
pert_seed <- c(57593, 19752, 09475, 39873, 62865, 48674, 44400, 27525, 59275, 82735, 20981)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVecList <- lapply(X = pert_real, FUN = wdsp.param, wNum = 4)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_9.csv")

# [10] Modify mean wind speed
pert_seed <- c(35903, 58267, 38956, 38665, 10993, 70862, 28622, 82695, 58926, 42394, 10989)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVecList <- lapply(X = pert_real, FUN = wdsp.param, wNum = 5)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_10.csv")

# [11] Modify SD wind speed
pert_seed <- c(83656, 97356, 07194, 65222, 88262, 11741, 49634, 69356, 71053, 92727, 29742)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVecList <- lapply(X = pert_real, FUN = wdsp.param, wNum = 6)
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_11.csv")

# [12] Modify max thistle density
# Note: reset misc_tDens before running the next code chunk
pert_seed <- c(39722, 36363, 01524, 39763, 01986, 69365, 57274, 78252, 32653, 60173, 10756)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
mVecList <- misc_tDens*pert_real
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_12.csv")

# [13] Modify ant nest density
# Note: reset misc_nDens before running the next code chunk
pert_seed <- c(98187, 92742, 82652, 51658, 01242, 39826, 49652, 27753, 60735, 35355, 49633)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
mVecList <- misc_nDens*pert_real
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_13.csv")

# [14] Modify max detection range from ant nests
# Note: reset misc_nRange before running the next code chunk
pert_seed <- c(29864, 10976, 69361, 09175, 81625, 50198, 34444, 29756, 79262, 39164, 90014)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
mVecList <- misc_nRange*pert_real
wv_pert <- c()
median(diff(wv_front))
write.csv(wv_pert, "Data/wv_pert_14.csv")

# Load in manually-combined version of CSVs above for future use
wv_pert <- read.csv("Data/wv_pert.csv")

# Perturbation notes below, replace 1 with rep number
# Use aVec <- aVecList[[1]] for parameters #1-5
# Use wVec <- wVecList[[1]] for parameters #6-11
# Use misc_tDens <- mVecList[[1]] for parameter #12
# Use misc_nDens <- mVecList[[1]] for parameter #13
# Use misc_nRange <- mVecList[[1]] for parameter #14





##### Run simulations for sensitivity analysis ------------------------------------------------------------

# Load simulation settings before running any wavespeed code
sets_plot <- FALSE     # (De)Activate wave plotting
sets_trim <- TRUE      # (De)Activate core area trimming
sets_years <- 1000     # Number of years to simulate
sets_trimD <- 1000     # Distance (m) behind wavefront to trim
sets_dCnts <- FALSE    # Track wind/ant dispersal contributions for each seed

# Set demography/dispersal parameters not covered in ADSP or WDSP
misc_tDens <- 10       # Max thistle density per metre
misc_nDens <- 0.2      # Ant nest density (nests/m)
misc_nRange <- 25      # Max detection range (m) from ant nests

# For each code chunk below, first run the lines from seed set to wVec definition (inclusive)
# Then run the entirety of 03_Wavespeeds and calculate median wavespeed
# Optional: save wavespeeds to CSV so there's no need to re-run code

# Baseline run
set.seed(96425)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_0.csv")

# [1] Modify head height intercept
set.seed(59735)
aVec <- adsp.param(aNum = 11, aVal = 0.8)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_1.csv")

# [2] Modify head height size-slope
set.seed(93462)
aVec <- adsp.param(aNum = 12, aVal = 0.8)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_2.csv")

# [3] Modify head height SD
set.seed(48165)
aVec <- adsp.param(aNum = 13, aVal = 0.8)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_3.csv")

# [4] Modify prob. of seed removal by ants
set.seed(81985)
aVec <- adsp.param(aNum = 20, aVal = 0.8)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_4.csv")

# [5] Modify prob. of predation if removed by ants
set.seed(20056)
aVec <- adsp.param(aNum = 21, aVal = 0.8)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_5.csv")

# [6] Modify prob. of seed release from capitulum
set.seed(79156)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 0.8)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_6.csv")

# [7] Modify surrounding vegetation height
set.seed(07195)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 2, wVal = 0.8)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_7.csv")

# [8] Modify mean seed terminal velocity
set.seed(59012)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 3, wVal = 0.8)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_8.csv")

# [9] Modify SD seed terminal velocity
set.seed(19754)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 4, wVal = 0.8)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_9.csv")

# [10] Modify mean wind speed
set.seed(49657)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 5, wVal = 0.8)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_10.csv")

# [11] Modify SD wind speed
set.seed(30567)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 6, wVal = 0.8)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_11.csv")

# [12] Modify max thistle density
# Note: reset misc_tDens before running any other code blocks
set.seed(41652)
misc_tDens <- misc_tDens*0.8
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_12.csv")

# [13] Modify ant nest density
# Note: reset misc_nDens before running any other code blocks
set.seed(63758)
misc_nDens <- misc_nDens*0.8
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_13.csv")

# [14] Modify max detection range from ant nests
# Note: reset misc_nRange before running any other code blocks
set.seed(90175)
misc_nRange <- misc_nRange*0.8
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)
median(diff(wv_front))
write.csv(wv_front, "Data/wv_elas_14.csv")

# Load in manually-combined version of CSVs above for future use
wv_elas <- read.csv("Data/wv_elas.csv")

