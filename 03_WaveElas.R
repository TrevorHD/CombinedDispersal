##### Set up functions for 1D expansion -------------------------------------------------------------------

# Function to simulate invasion wave
wave.sim <- function(dVec, sVec){
  
  # Set up framework for parallel processing
  cl <- makeCluster(detectCores() - 1, type = "PSOCK")
  parallel::setDefaultCluster(cl)
  clusterEvalQ(cl, {library(SuppDists)
                    library(MASS)
                    library(tidyverse)
                    library(truncnorm)})
    
  # Function to see if a seed is taken to the nearest nest
  nestsearch <- function(d, range, sVec){
    dists <- abs(d - nestsR)
    centre <- nestsR[which.min(dists)]
    toProb <- ant(min(dists), sVec)
    toNest <- sample(c(0, 1), 1, prob = c(1 - toProb, toProb))
    ifelse(toNest == 1 && min(dists) <= range, return(centre), return(d))}
  
  # Estimate dispersal distances from given point
  kern <- function(n, h, species, sVec, d0 = 0){
    d <- wald(n, h, species, sVec) + d0
    return(d)}
  
  #disp_index <- rep(1, length(opt_ht))
  
  # NEW: Estimate dispersal distances from given point
  #kern <- function(n, h, species, sVec, disp = data_disp, d0 = 0){
    
    # Pull column index based on height
    #ci <- h*100 - 15
    
    # Pull row index for given column
    #index_h <- disp_index[ci]
    
    # Pull index:(index + n) from column and save distances, then advance index by n
    # Also accounts for wrap-around
    #if(index_h + n > 1000000){
    #  dists <- c(disp[index_h:1000000, ci], disp[0:(index_h + n - 1000001), ci]) + d0
    #  disp_index[ci] <<- index_h + n - 1000000
    #} else {
    #  dists <- disp[index_h:(index_h + n - 1), ci] + d0
    #  disp_index[ci] <<- index_h + n}
    
    # Return distances
    #return(dists)}
  
  # Choose species to model
  species <- "CN"
  
  # Set various parameters for wave model
  nestOn <- TRUE    # Should ant nests be included
  range <- 15       # Max detection range (m) from ant nests
  nDens <- 0.2      # Ant nest density (nests/m)
  nYear <- 1000     # Number of years to simulate
  trim <- TRUE      # Should core area of wave be trimmed?
  trimAmt <- 1000   # Distance (m) behind wavefront to trim
  tDens <- 5        # Max thistle density per metre
  plotOn <- FALSE   # Plot wave?
  
  # Run invasion wave simulation
  for(i in 1:nYear){
    
    # Initialise simulation data
    if(i == 1){
      
      # Generate nests
      nests <- sample(seq(0, nYear*200, by = 0.1), nDens*nYear*200) + 0.01
      
      # Initialise data; start with a single rosette
      plants <- data.frame(d = 0.01, stage = 0, rsize = demo("rsize", species, dVec, n = 1))
      seedsAG <- data.frame(matrix(ncol = 2, nrow = 0))
      colnames(seedsAG) <- c("d", "germ")
      seedsSB <- data.frame(matrix(ncol = 2, nrow = 0))
      colnames(seedsSB) <- c("d", "germ")
      vals <- c()}
    
    # Surviving above-ground seeds become rosettes
    seedsAG$stage <- rep(0, nrow(seedsAG))
    seedsAG$rsize <- demo("rsize", species, dVec, n = nrow(seedsAG))
    seedsAG <- seedsAG[, !names(seedsAG) == c("germ")]
    plants <- rbind(plants, seedsAG)
    
    # Kill rosettes if density is too high
    # Do this by sorting so that adults come first and are prioritised
    if(nrow(plants) > 0){
      plants %>% 
        mutate(bin = as.integer(cut(d, breaks = seq(0, ceiling(max(plants$d)), 1)))) %>% 
        arrange(bin, desc(stage)) %>% 
        group_by(bin) %>%
        slice(1:pmin(tDens, n())) %>% 
        data.frame() -> plants
      plants <- plants[, !names(plants) == c("bin")]}
    
    # Reset seeds (except for seed bank)
    seedsAG <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsAG) <- c("d", "germ")
    seedsNew <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsNew) <- c("d", "germ")
    
    # Plot density over space
    if(i == 1){
      PlotList <- list()}
    if(plotOn == TRUE){
      PlotList[[i]] <- plants$d}
    
    # Simulate rosette survival and remove rosettes that do not survive
    # Never kill when only one rosette; prevents code from crashing
    if(nrow(plants) > 1){
      plants <- plants[demo("survival", species, dVec, rsize = plants$rsize) == 1, ]}
    
    # Trim core areas as wave progresses to save computational resources
    if(trim == TRUE & nrow(plants) > 0){
      plants <- plants[plants$d > max(plants$d) - trimAmt, ]
      seedsSB <- seedsSB[seedsSB$d > max(plants$d) - trimAmt, ]
      nestsR <- nests[nests > min(plants$d) & nests < max(plants$d) + 1000]}
    
    # Simulate bolting and flowering; individuals that don't will remain rosettes
    plants$stage <- demo("flowering", species, dVec, rsize = plants$rsize)
    
    # Survival and establishment of seeds already in seed bank
    if(nrow(seedsSB) > 0){
      seedsSB <- seedsSB[demo("surSB", species, dVec, n = nrow(seedsSB)) == 1, ]
      if(nrow(seedsSB) > 0){
        vec <- demo("estSB", species, dVec, n = nrow(seedsSB))
        temp1 <- seedsSB[vec == 1, ]
        temp1$stage <- rep(0, nrow(temp1))
        temp1$rsize <- demo("rsize", species, dVec, n = nrow(temp1))
        temp1 <- temp1[, !names(temp1) == c("germ")]}
      plants <- na.omit(rbind(plants, temp1))
      seedsSB <- seedsSB[vec == 0, ]}
    
    # Simulate dispersal from flowering adults
    if(sum(plants$stage == 1) > 0){
      temp2 <- plants[plants$stage == 1, ]
      nflow <- demo("flowers", species, dVec, rsize = temp2$rsize)
      f1 <- unlist(sapply(nflow, demo, dType = "height", species = species, dVec = dVec))
      s1 <- rep(demo("seeds", "CN", dVec), times = sum(nflow))
      d1 <- unlist(mapply(x = temp2$d, times = nflow, rep))
      clusterExport(cl, c("kern", "wald", "s1", "f1", "d1", "species", "sVec"))
      seedsNew <- unlist(clusterMap(cl, kern, n = s1, h = f1, d0 = d1,
                                    MoreArgs = list(species = species, sVec = sVec)))
      seedsNew <- data.frame(cbind(seedsNew, rep(0, length(seedsNew))))
      names(seedsNew) <- c("d", "germ")}
    
    # Simulate aboveground establishment from dispersed seeds
    seedsNew$germ <- demo("estAG", species, dVec, n = nrow(seedsNew))
    seedsAG <- na.omit(rbind(seedsAG, seedsNew[seedsNew$germ == 1, ]))
    seedsNew <- seedsNew[seedsNew$germ == 0, ]
    
    # Entry of non-establishing seeds into seed bank
    if(nrow(seedsNew) > 0){
      seedsNew$germ <- demo("entSB", species, dVec, n = nrow(seedsNew))
      seedsSB <- na.omit(rbind(seedsSB, seedsNew[seedsNew$germ == 1, ]))
      seedsSB$germ <- rep(0, nrow(seedsSB))}
    
    # Simulate secondary seed dispersal via ants
    if(nestOn == TRUE){
      clusterExport(cl, c("seedsAG", "nestsearch", "ant", "nestsR", "range", "sVec"))
      seedsAG$d <- parSapply(cl, seedsAG$d, nestsearch, range = range, sVec = sVec)}
    
    # Store wavefront distance
    vals <- c(vals, max(plants$d))
    
    # Kill all adults after they reproduce
    plants <- plants[plants$stage != 1, ]
    
    # Print procedure progress
    shell("cls")
    cat(paste0(i, "/", nYear, " (", round(i/nYear, 3)*100, "%) complete"))}
  
  # Return list of wavefront positions and all plant positions
  return(list(wavefront = vals, positions = PlotList))}

# Calculate baseline wavespeeds
wave.base <- function(){
  
  # Start system timer
  time.start <- Sys.time()
  
  # Initialise parameter vectors
  dVec1 <- demo.param(dNum = 1, dVal = 1)
  sVec1 <- wald.param(sNum = 1, sVal = 1)
  
  # Run simulation without any parameter changes
  wave1 <- wave.sim(dVec1, sVec1)
  
  # Calculate mean wavespeed
  mWave1 <- mean(diff(wave1$wavefront))
  
  # Calculate procedure time
  time.elapsed <- as.numeric(difftime(Sys.time(), time.start, units = "hours"))
  
  # Return list of calculated quantities
  return(list(procedure.time = paste0(time.elapsed, " hours"), wSpeedMean1 = mWave1,
              wSpeed1 = diff(wave1$wavefront), wFront1 = wave1$wavefront))}

# Calculate wavespeed elasticity
wave.elas <- function(dNum, dVal, sNum, sVal, baseObj){
  
  # Let dNum be the demographic parameter number
  # Let dVal be the proportion to multiply the original parameter by
  # Let baseObj be the baseline wave object
  
  # Start system timer
  time.start <- Sys.time()
  
  # Initialise parameter vectors
  if(dVal == 1 & sVal != 1){
    dVec2 <- demo.param(dNum = 1, dVal = 1)
    sVec2 <- wald.param(sNum = sNum, sVal = sVal)}
  if(dVal != 1 & sVal == 1){
    dVec2 <- demo.param(dNum = dNum, dVal = dVal)
    sVec2 <- wald.param(sNum = 1, sVal = 1)}
  
  # Run simulation with increase/decrease on a specified parameter
  wave2 <- wave.sim(dVec2, sVec2)
  
  # Calculate mean wavespeeds
  mWave1 <- baseObj$wSpeedMean1
  mWave2 <- mean(diff(wave2$wavefront))
  
  # Calculate percent change in wavespeed and parameter
  wPct <- (mWave2 - mWave1)/mWave1*100
  if(dVal == 1 & sVal != 1){
    pPct <- (sVal - 1)*100}
  if(dVal != 1 & sVal == 1){
    pPct <- (dVal - 1)*100}
  
  # Calculate elasticity
  elas <- wPct/pPct
  
  # Calculate procedure time
  time.elapsed <- as.numeric(difftime(Sys.time(), time.start, units = "hours"))
  
  # Return list of calculated quantities
  return(list(procedure.time = paste0(time.elapsed, " hours"), elasticity = elas,
              wSpeedMean1 = mWave1, wSpeedMean2 = mWave2,
              wSpeed1 = baseObj$wSpeed1, wSpeed2 = diff(wave2$wavefront),
              wFront1 = baseObj$wFront1, wFront2 = wave2$wavefront))}

# Generate GIF of population spread
# Limit to 10000 m (nYear = 100 recommended)
generatePlots <- function(type = "hist"){
  for(i in 1:length(wave1[[2]])){
    lower <- c()
    positions <- unlist(wave1[[2]][i])
    if(length(positions) == 0){
      positions <- 0}
    if(max(positions > 500)){
      minBin <- sum(positions > floor(min(positions)) & positions < ceiling(min(positions)))
      lower <- c(rep(0:(floor(min(positions) - 1)), 10) + 0.01,
                 rep(floor(min(positions)) + 0.01, 10 - minBin))
      if(min(positions) < 1){
        lower <- sort(lower)[-c(1:20)]}}
    if(type == "hist"){
      hist(c(lower, positions), breaks = seq(0, 10000, by = 20), xlim = c(0, 10000), ylim = c(0, 250),
           xaxt = "n", yaxt = "n", xlab = "Distance (m)", ylab = "Count", main = "")
      axis(1, at = seq(0, 10000, 2000))
      axis(2, at = seq(0, 250, 50))
      text(x = 9400, y = 230, paste0("t = ", i))
      box()}
    if(type == "density"){
      plotdata <- hist(c(lower, positions), breaks = seq(0, 10000, by = 20), plot = FALSE)
      plot(plotdata$mids, plotdata$density/max(plotdata$density),
           type = "l", xlim = c(0, 10000), ylim = c(0, 1.25),
           xaxt = "n", yaxt = "n", xlab = "Distance (m)", ylab = "Relative Density", main = "")
      axis(1, at = seq(0, 10000, 2000))
      axis(2, at = seq(0, 1.25, 0.25))
      text(x = 9400, y = 1.15, paste0("t = ", i))}}}
save_gif(generatePlots("hist"), "Spread1.gif", delay = 0.3, width = 1280, height = 720, res = 144)
save_gif(generatePlots("density"), "Spread2.gif", delay = 0.3, width = 1280, height = 720, res = 144)





##### Simulate wavespeeds ---------------------------------------------------------------------------------

# Set seed to ensure replicability
set.seed(8675309)

# Set counter for dispersal height index
disp_index <- rep(1, length(opt_ht))

# Test run code
wave.base()
