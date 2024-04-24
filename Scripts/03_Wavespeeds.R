# Start system timer
wv_time <- Sys.time()

# Set up cluster for parallel processing
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
parallel::setDefaultCluster(cl)
clusterEvalQ(cl, {library(SuppDists)
  library(MASS)
  library(tidyverse)
  library(truncnorm)})

# Enforce current RNG seed for all parallel processes
clusterSetRNGStream(cl = cl)

# Run invasion wave simulation
for(i in 1:(sets_years + 1)){
  
  # Initialise simulation data
  if(i == 1){
    
    # Generate nests
    nests <- sample(seq(0, sets_years*200, by = 0.1), misc_nDens*sets_years*200) + 0.01
    nestsR <- nests
    
    # Initialise data; start with a single rosette and no seeds
    plants <- data.frame(d = 0.01, stage = 0, rsize = adsp.demo("size", aVec, n = 1))
    seedsNew <- data.frame(matrix(ncol = 2, nrow = 0))
    seedsEst <- data.frame(matrix(ncol = 2, nrow = 0))
    seedsSB <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsNew) <- c("d", "germ")
    colnames(seedsEst) <- c("d", "germ")
    colnames(seedsSB) <- c("d", "germ")
    wv_front <- c()}
  
  # Simulate rosette survival and remove rosettes that do not survive
  if(nrow(plants) > 0){
    plants <- plants[adsp.demo("surv", aVec, n = nrow(plants)) == 1, ]}
  
  # Simulate rosette growth
  if(nrow(plants) > 0){
    plants$rsize = adsp.demo("grow", aVec, rsize = plants$rsize)}
  
  # Simulate bolting and flowering; individuals that don't will remain rosettes
  if(nrow(plants) > 0){
    plants$stage <- adsp.demo("flow", aVec, n = nrow(plants))}
  
  # Estimate number of flower heads, then get the number of viable seeds and their release heights
  # Viable means surviving pre-/post-dispersal predation, and establishing or entering seed bank
  if(sum(plants$stage == 1) > 0){
    temp1 <- plants[plants$stage == 1, ]
    nflow <- adsp.demo("head", aVec, rsize = temp1$rsize)
    f1 <- adsp.demo("hdht", aVec, rsize = rep(temp1$rsize, times = nflow))
    s1 <- rep(adsp.demo("seed", aVec), times = sum(nflow))
    d1 <- rep(temp1$d, times = nflow)}
  
  # Simulate primary dispersal of viable seeds via wind
  if(sum(s1) > 0){
    clusterExport(cl, c("wdsp.disp", "wdsp.wald", "s1", "f1", "d1", "wVec"))
    seedsNew <- unlist(clusterMap(cl, wdsp.disp, n = s1, H = f1, d0 = d1,
                                  MoreArgs = list(wVec = wVec)))
    seedsNew <- data.frame(cbind(seedsNew, rep(0, length(seedsNew))))
    names(seedsNew) <- c("d", "germ")
    seedsNew <- seedsNew[seedsNew$d >= 0, ]}
  
  # Simulate secondary seed dispersal via ants
  if(nrow(seedsNew) > 0 & sets_nest == TRUE){
    temp2 <- adsp.demo("ants", aVec, n = nrow(seedsNew))
    temp3 <- seedsNew[temp2 == 1, ]
    temp4 <- seedsNew[temp2 == 0, ]
    clusterExport(cl, c("adsp.disp", "nestsR", "misc_nRange", "temp3"))
    temp3$d <- unlist(parSapply(cl, temp3$d, adsp.disp, range = misc_nRange))
    seedsNew <- rbind(temp3, temp4)}
  
  # Simulate establishment from dispersed seeds not in seed bank
  if(nrow(seedsNew) > 0){
    seedsNew$germ <- adsp.demo("estb", aVec, n = nrow(seedsNew))
    seedsEst <- na.omit(rbind(seedsEst, seedsNew[seedsNew$germ == 1, ]))
    seedsNew <- seedsNew[seedsNew$germ == 0, ]}
  
  # Simulate survival, then establishment, of seeds already in seed bank
  if(nrow(seedsSB) > 0){
    seedsSB <- seedsSB[adsp.demo("SBsurv", aVec, n = nrow(seedsSB)) == 1, ]
    if(nrow(seedsSB) > 0){
      vec <- adsp.demo("SBestb", aVec, n = nrow(seedsSB))
      temp5 <- seedsSB[vec == 1, ]
      temp5$germ <- rep(1, nrow(temp5))}
    if(nrow(seedsSB) > 0){
      seedsSB <- seedsSB[vec == 0, ]
      seedsEst <- na.omit(rbind(seedsEst, temp5))}}
  
  # All remaining seeds not experiencing post-predation death enter seed bank
  # Then reset all seeds already assigned to establishment or seed bank
  if(nrow(seedsNew) > 0){
    seedsNew$germ <- rep(0, nrow(seedsNew))
    seedsSB <- na.omit(rbind(seedsSB, seedsNew))}
  seedsNew <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(seedsNew) <- c("d", "germ")
  
  # Kill adults that have bolted and/or reproduced
  plants <- plants[plants$stage == 0, ]
  
  # Simulate rosette size of established individuals
  # Then reset seeds that have been assigned to establishment
  if(nrow(seedsEst) > 0){
    seedsEst$stage <- rep(0, nrow(seedsEst))
    seedsEst$rsize <- adsp.demo("size", aVec, n = nrow(seedsEst))
    seedsEst <- seedsEst[, !names(seedsEst) == c("germ")]
    plants <- rbind(plants, seedsEst)
    seedsEst <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsEst) <- c("d", "germ")}
  
  # Kill any rosettes in excess of maximum density per metre
  # Larger rosettes get priority when selecting which ones are kept
  if(misc_tDens == 0){
    plants <- plants[-c(1:nrow(plants)), ]}
  if(nrow(plants) > 0){
    plants %>% 
      filter(d > 0) %>%
      mutate(bin = as.integer(cut(d, breaks = seq(0, ceiling(max(plants$d)), 1)), right = FALSE)) %>% 
      arrange(bin, d) %>% 
      group_by(bin) %>%
      slice(1:pmin(misc_tDens, n())) %>% 
      data.frame() -> plants
    plants <- plants[, !names(plants) == c("bin")]}
  
  # If no plants remain (unlikely), reset at origin so that simulation does not fail
  if(nrow(plants) < 1){
    plants <- data.frame(d = 0.01, stage = 0, rsize = adsp.demo("size", aVec, n = 1))}
  
  # Store wavefront distance
  wv_front <- c(wv_front, max(plants$d))
  
  # Plot density over space
  if(i == 1){
    wv_plots <- list()}
  if(sets_plot == TRUE){
    wv_plots[[i]] <- plants$d}
  
  # Trim core areas as wave progresses to save computational resources
  if(sets_trim == TRUE & nrow(plants) > 0){
    plants <- plants[plants$d > max(plants$d) - sets_trimD, ]
    seedsSB <- seedsSB[seedsSB$d > max(plants$d) - sets_trimD, ]
    nestsR <- nests[nests > min(plants$d) & nests < max(plants$d) + sets_trimD]}
  
  # Print procedure progress
  shell("cls")
  cat(paste0(i - 1, "/", sets_years, " (", round(i/(sets_years + 1), 3)*100, "%) complete with ",
             round(as.numeric(difftime(Sys.time(), wv_time, units = "mins")), 2),
             " minutes elapsed"))}

# Close cluster
stopCluster(cl)

# Remove variables no longer in use
remove(plants, seedsEst, seedsNew, seedsSB, temp1, temp2, temp3, temp4, temp5,
       d1, f1, i, nests, nestsR, nflow, s1, vec)

# Print final completion message
shell("cls")
cat(paste0("Procedure complete with ",
           round(as.numeric(difftime(Sys.time(), wv_time, units = "mins")), 2),
           " minutes elapsed"))

