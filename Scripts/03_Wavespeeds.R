##### Simulate 1D expansion -------------------------------------------------------------------------------

# Set changes to base parameters (1 if not doing sensitivity analysis)
dVec <- demo.param(dNum = 1, dVal = 1)
sVec <- wald.param(sNum = 1, sVal = 1)

# Start system timer
wv_time <- Sys.time()

# Set up cluster for parallel processing
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
parallel::setDefaultCluster(cl)
clusterEvalQ(cl, {library(SuppDists)
  library(MASS)
  library(tidyverse)
  library(truncnorm)})

# Set various parameters for wave model
nestOn <- TRUE    # Should ant nests be included
range <- 15       # Max detection range (m) from ant nests
nDens <- 0.2      # Ant nest density (nests/m)
nYear <- 1000     # Number of years to simulate
trim <- TRUE      # Should core area of wave be trimmed?
trimAmt <- 1000   # Distance (m) behind wavefront to trim
tDens <- 10       # Max thistle density per metre
plotOn <- TRUE    # Plot wave?

# Run invasion wave simulation
for(i in 1:nYear){
  
  # Initialise simulation data
  if(i == 1){
    
    # Generate nests
    nests <- sample(seq(0, nYear*200, by = 0.1), nDens*nYear*200) + 0.01
    
    # Initialise data; start with a single rosette
    plants <- data.frame(d = 0.01, stage = 0, rsize = demo("rsize", dVec, n = 1))
    seedsAG <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsAG) <- c("d", "germ")
    seedsSB <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsSB) <- c("d", "germ")
    wv_front <- c()}
  
  # Simulate surviving above-ground seeds becoming rosettes
  seedsAG$stage <- rep(0, nrow(seedsAG))
  seedsAG$rsize <- demo("rsize", dVec, n = nrow(seedsAG))
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
    wv_plots <- list()}
  if(plotOn == TRUE){
    wv_plots[[i]] <- plants$d}
  
  # Simulate rosette survival and remove rosettes that do not survive
  # Never kill when only one plant, and reset another rosette at origin if no plants remain
  # Both of these measures prevent the simulation from failing unexpectedly
  if(nrow(plants) > 1){
    plants <- plants[demo("survival", dVec, rsize = plants$rsize) == 1, ]}
  if(nrow(plants) == 0)
    plants <- data.frame(d = 0.01, stage = 0, rsize = demo("rsize", dVec, n = 1))
  
  # Simulate rosette growth
  if(nrow(plants) > 0){
    plants$rsize = demo("grow", dVec, rsize = plants$rsize)}
  
  # Trim core areas as wave progresses to save computational resources
  if(trim == TRUE & nrow(plants) > 0){
    plants <- plants[plants$d > max(plants$d) - trimAmt, ]
    seedsSB <- seedsSB[seedsSB$d > max(plants$d) - trimAmt, ]
    nestsR <- nests[nests > min(plants$d) & nests < max(plants$d) + 1000]}
  
  # Simulate bolting and flowering; individuals that don't will remain rosettes
  plants$stage <- demo("flowering", dVec, rsize = plants$rsize)
  
  # Simulate survival, then establishment, of seeds already in seed bank
  if(nrow(seedsSB) > 0){
    seedsSB <- seedsSB[demo("surSB", dVec, n = nrow(seedsSB)) == 1, ]
    if(nrow(seedsSB) > 0){
      vec <- demo("estSB", dVec, n = nrow(seedsSB))
      temp1 <- seedsSB[vec == 1, ]
      temp1$stage <- rep(0, nrow(temp1))
      temp1$rsize <- demo("rsize", dVec, n = nrow(temp1))
      temp1 <- temp1[, !names(temp1) == c("germ")]}
    if(nrow(seedsSB) > 0){
      seedsSB <- seedsSB[vec == 0, ]
      plants <- na.omit(rbind(plants, temp1))}}
  
  # Simulate dispersal from flowering adults
  if(sum(plants$stage == 1) > 0){
    temp2 <- plants[plants$stage == 1, ]
    nflow <- demo("flowers", dVec, rsize = temp2$rsize)
    f1 <- unlist(sapply(nflow, demo, dType = "height", dVec = dVec))
    s1 <- rep(demo("seeds", dVec), times = sum(nflow))
    d1 <- unlist(mapply(x = temp2$d, times = nflow, rep))
    clusterExport(cl, c("kern", "wald", "s1", "f1", "d1", "sVec"))
    seedsNew <- unlist(clusterMap(cl, kern, n = s1, h = f1, d0 = d1,
                                  MoreArgs = list(sVec = sVec)))
    seedsNew <- data.frame(cbind(seedsNew, rep(0, length(seedsNew))))
    names(seedsNew) <- c("d", "germ")
    seedsNew <- seedsNew[seedsNew$d >= 0, ]}
  
  # Simulate aboveground establishment from dispersed seeds
  seedsNew$germ <- demo("estAG", dVec, n = nrow(seedsNew))
  seedsAG <- na.omit(rbind(seedsAG, seedsNew[seedsNew$germ == 1, ]))
  seedsNew <- seedsNew[seedsNew$germ == 0, ]
  
  # Simulate entry of non-establishing seeds into seed bank
  if(nrow(seedsNew) > 0){
    seedsNew$germ <- demo("entSB", dVec, n = nrow(seedsNew))
    seedsSB <- na.omit(rbind(seedsSB, seedsNew[seedsNew$germ == 1, ]))
    seedsSB$germ <- rep(0, nrow(seedsSB))}
  
  # Simulate secondary seed dispersal via ants
  if(nestOn == TRUE){
    clusterExport(cl, c("seedsAG", "nestsearch", "ant", "nestsR", "range", "sVec"))
    seedsAG$d <- parSapply(cl, seedsAG$d, nestsearch, range = range, sVec = sVec)}
  
  # Store wavefront distance
  wv_front <- c(wv_front, max(plants$d))
  
  # Kill all adults after they reproduce
  plants <- plants[plants$stage != 1, ]
  
  # Print procedure progress
  shell("cls")
  cat(paste0(i, "/", nYear, " (", round(i/nYear, 3)*100, "%) complete with ",
             round(as.numeric(difftime(Sys.time(), wv_time, units = "mins")), 2),
             " minutes elapsed"))}

# Close cluster
stopCluster(cl)

# Remove variables no longer in use
remove(plants, seedsAG, seedsNew, seedsSB, temp1, temp2, d1, f1, i, nDens, nestOn,
       nests, nestsR, nflow, nYear, plotOn, range, s1, tDens, trim, trimAmt, vec)

# Print final completion message
shell("cls")
cat(paste0("Procedure complete with ",
           round(as.numeric(difftime(Sys.time(), wv_time, units = "mins")), 2),
           " minutes elapsed"))

