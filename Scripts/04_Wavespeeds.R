##### Run population spread model -------------------------------------------------------------------------

# Start system timer
wv_time <- Sys.time()

# Set up cluster for parallel processing
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
parallel::setDefaultCluster(cl)
clusterEvalQ(cl, {library(SuppDists); library(MASS);
                  library(tidyverse); library(truncnorm)})

# Enforce current RNG seed for all parallel processes
clusterSetRNGStream(cl = cl)

# Run wave simulation for specified number of years
for(i in 1:(sets_years + 1)){
  
  # Initialise simulation data
  if(i == 1){
    
    # Generate ant nests at the specified average density
    if(misc_nDens > 0){
      nests <- sample(seq(0, 1/misc_nDens - 0.01, by = 0.01), misc_nDens*sets_years*200, replace = TRUE)
      nests <- round(nests + seq(0, sets_years*200 - 1/misc_nDens, by = 1/misc_nDens), 2)}
    if(misc_nDens == 0){
      nests <- c()}
    
    # Generate vector to store wavefront positions
    wv_front <- c()
    
    # Initialise data, starting with a single rosette (random size) and no seeds
    # If tracking dispersal details, initialise same data as above, but with an ID column
    if(sets_dCnts == FALSE){
      plants <- data.frame(d = 0.01, stage = 0, rsize = adsp.demo("size", aVec, n = 1))
      seedsNew <- data.frame(matrix(ncol = 2, nrow = 0))
      seedsEst <- data.frame(matrix(ncol = 2, nrow = 0))
      seedsSB <- data.frame(matrix(ncol = 2, nrow = 0))
      colnames(seedsNew) <- c("d", "germ")
      colnames(seedsEst) <- c("d", "germ")
      colnames(seedsSB) <- c("d", "germ")}
    if(sets_dCnts == TRUE){
      plants <- data.frame(ID = 0, d = 0.01, stage = 0, rsize = adsp.demo("size", aVec, n = 1))
      seedsNew <- data.frame(matrix(ncol = 3, nrow = 0))
      seedsEst <- data.frame(matrix(ncol = 3, nrow = 0))
      seedsSB <- data.frame(matrix(ncol = 3, nrow = 0))
      colnames(seedsNew) <- c("ID", "d", "germ")
      colnames(seedsEst) <- c("ID", "d", "germ")
      colnames(seedsSB) <- c("ID", "d", "germ")}
  
    # If tracking dispersal details, generate data frames to store various dispersal stats
    if(sets_dCnts == TRUE){
      wv_ids <- data.frame(matrix(ncol = 3, nrow = 0))
      wv_dCnts <- data.frame(matrix(ncol = 2, nrow = 0))
      wv_disp_w <- data.frame(matrix(ncol = 7, nrow = 0))
      wv_disp_a <- data.frame(matrix(ncol = 7, nrow = 0))
      colnames(wv_ids) <- c("ID", "wd", "ad")
      colnames(wv_dCnts) <- c("Wind", "Ant")
      colnames(wv_disp_w) <- c("MeanAbs", "MedianAbs", "MaxAbs", "25P", "10P", "5P", "1P")
      colnames(wv_disp_a) <- c("MeanAbs", "MedianAbs", "MaxAbs", "25P", "10P", "5P", "1P")}}
    
  # Simulate rosette survival and remove rosettes that do not survive
  # To prevent unexpected simulation failure, do this only if number of individuals exceeds 1
  if(nrow(plants) > 1){
    plants <- plants[adsp.demo("surv", aVec, n = nrow(plants)) == 1, ]}
  
  # Simulate rosette growth
  if(nrow(plants) > 0){
    plants$rsize = adsp.demo("grow", aVec, rsize = plants$rsize)}
  
  # Simulate bolting and flowering; individuals that don't will remain rosettes
  if(nrow(plants) > 0){
    plants$stage <- adsp.demo("flow", aVec, n = nrow(plants))}
  
  # Estimate number of flower heads, then get the number of seeds and their release heights
  # Keep only seeds that survive pre-/post-dispersal predation, and establish or enter seed bank
  if(sum(plants$stage == 1) > 0){
    temp1 <- plants[plants$stage == 1, ]
    nflow <- adsp.demo("head", aVec, rsize = temp1$rsize)
    f1 <- adsp.demo("hdht", aVec, rsize = rep(temp1$rsize, times = nflow))
    s1 <- rep(adsp.demo("seed", aVec), times = sum(nflow))
    d1 <- rep(temp1$d, times = nflow)
    remove(temp1, nflow)}
  if(sum(plants$stage == 1) == 0){
    s1 <- 0}
  
  # Simulate primary dispersal of viable seeds via wind, recording new seed locations
  # If tracking dispersal details, record wind dispersal distance for each seed
  if(sum(s1) > 0){
    clusterExport(cl, c("wdsp.disp", "wdsp.wald", "s1", "f1", "d1", "wVec"))
    seedsNew <- unlist(clusterMap(cl, wdsp.disp, n = s1, H = f1, d0 = d1, MoreArgs = list(wVec = wVec)))
    seedsNew <- data.frame(d = round(seedsNew, 2), germ = rep(0, length(seedsNew)))
    if(sets_dCnts == TRUE){
      temp2 <- seedsNew$d - rep(d1, times = s1)}
    remove(f1, d1)}
  
  # Remove seeds that fall below minimum boundary of simulation area
  if(sum(s1) > 0){
    if(sets_dCnts == TRUE){
      temp3 <- temp2[seedsNew$d >= min(plants$d)]}
    seedsNew <- seedsNew[seedsNew$d >= min(plants$d), ]}
  
  # If tracking dispersal details, add IDs and match wind dispersal distances to each seed
  # Then calculate select statistics for seeds experiencing wind dispersal (forward only)
  if(sum(s1) > 0){
    if(sets_dCnts == TRUE){
      seedsNew$ID <- wdsp.sIDs(plants$ID, seedsSB$ID, nrow(seedsNew))
      wv_ids_wd <- data.frame(ID = seedsNew$ID, wd = temp3)
      temp2 <- temp2[temp2 > 0]
      wv_disp_w <- rbind(wv_disp_w, c(mean(temp2), median(temp2), max(temp2), quantile(temp2, 0.75),
                                      quantile(temp2, 0.90), quantile(temp2, 0.95), quantile(temp2, 0.99)))
      names(wv_disp_w) <- c("MeanAbs", "MedianAbs", "MaxAbs", "25P", "10P", "5P", "1P")
      wv_disp_w <- round(wv_disp_w, 2)
      remove(temp2, temp3)}
    remove(s1)}
  
  # Limit ant search to current simulation area to save computational resources
  nestsR <- nests[nests > min(plants$d) & nests < max(seedsNew$d) + 2*1/misc_nDens]
  
  # Simulate secondary seed dispersal via ants, recording new seed locations
  # If tracking dispersal details, record ant dispersal distances for each seed 
  if(nrow(seedsNew) > 0){
    if(misc_nDens != 0 & misc_nRange != 0){
      temp4 <- adsp.demo("ants", aVec, n = nrow(seedsNew))
      temp5 <- seedsNew$d
      temp6 <- seedsNew[temp4 == 1, ]$d
      clusterExport(cl, c("adsp.disp", "nestsR", "misc_nRange", "temp6"))
      temp6 <- unlist(parSapply(cl, temp6, adsp.disp, range = misc_nRange, nestList = nestsR))
      seedsNew[temp4 == 1, ]$d <- temp6
      if(sets_dCnts == TRUE){
        temp7 <- seedsNew$d - temp5}
      remove(temp4, temp5, temp6)}
    if(misc_nDens == 0 | misc_nRange == 0){
      if(sets_dCnts == TRUE){
        temp7 <- rep(0, nrow(seedsNew))}}}

  # If tracking dispersal details, match ant dispersal distances to each seed
  # Then calculate select statistics for seeds experiencing ant dispersal (forward only)
  if(nrow(seedsNew) > 0){
    if(sets_dCnts == TRUE){
      wv_ids_ad <- data.frame(ID = seedsNew$ID, ad = temp7)
      wv_ids <- rbind(wv_ids, merge(wv_ids_wd, wv_ids_ad))
      if(misc_nDens != 0 & misc_nRange != 0){
        temp7 <- temp7[temp7 > 0]}
      wv_disp_a <- rbind(wv_disp_a, c(mean(temp7), median(temp7), max(temp7), quantile(temp7, 0.75),
                                      quantile(temp7, 0.90), quantile(temp7, 0.95), quantile(temp7, 0.99)))
      names(wv_disp_a) <- c("MeanAbs", "MedianAbs", "MaxAbs", "25P", "10P", "5P", "1P")
      wv_disp_a <- round(wv_disp_a, 2)
      remove(temp7)}}
  
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
      temp8 <- seedsSB[vec == 1, ]
      temp8$germ <- rep(1, nrow(temp8))}
    if(nrow(seedsSB) > 0){
      seedsSB <- seedsSB[vec == 0, ]
      seedsEst <- na.omit(rbind(seedsEst, temp8))
      remove(vec, temp8)}}
  
  # Simulate entry into seed bank for remaining seeds
  # Then reset "new seeds" after assignment to establishment or seed bank
  if(nrow(seedsNew) > 0){
    seedsNew$germ <- rep(0, nrow(seedsNew))
    seedsSB <- na.omit(rbind(seedsSB, seedsNew))}
  if(sets_dCnts == FALSE){
    seedsNew <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsNew) <- c("d", "germ")}
  if(sets_dCnts == TRUE){
    seedsNew <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(seedsNew) <- c("ID", "d", "germ")}
  
  # Remove adults that have bolted and/or reproduced
  plants <- plants[plants$stage == 0, ]
  
  # Simulate rosette size of established individuals
  # Move established individuals into plants data frame
  if(nrow(seedsEst) > 0){
    seedsEst$stage <- rep(0, nrow(seedsEst))
    seedsEst$rsize <- adsp.demo("size", aVec, n = nrow(seedsEst))
    seedsEst <- seedsEst[, !names(seedsEst) == c("germ")]
    plants <- rbind(plants, seedsEst)}
  if(sets_dCnts == FALSE){
    seedsEst <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(seedsEst) <- c("d", "germ")}
  if(sets_dCnts == TRUE){
    seedsEst <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(seedsEst) <- c("ID", "d", "germ")}
  
  # Remove any rosettes in excess of maximum density per metre
  # Larger rosettes get priority when selecting which ones are kept
  if(misc_tDens == 0){
    plants <- plants[-c(1:nrow(plants)), ]}
  if(nrow(plants) > 0){
    plants %>% 
      filter(d > 0) %>%
      mutate(bin = as.integer(cut(d, breaks = seq(0, ceiling(max(plants$d)), 1)), right = FALSE)) %>% 
      arrange(bin, desc(rsize), d) %>% 
      group_by(bin) %>%
      slice(1:pmin(misc_tDens, n())) %>% 
      data.frame() -> plants
    plants <- plants[, !names(plants) == c("bin")]}
  
  # Reset at origin if no plants remain (unlikely) so that simulation does not fail
  if(nrow(plants) < 1){
    plants <- data.frame(d = 0.01, stage = 0, rsize = adsp.demo("size", aVec, n = 1))}
  
  # Trim core areas as wave progresses to save computational resources
  if(nrow(plants) > 0){
    if(sets_trim == TRUE){
      plants <- plants[plants$d > max(plants$d) - sets_trimD, ]
      seedsSB <- seedsSB[seedsSB$d > max(plants$d) - sets_trimD, ]}}
  
  # Store wavefront distance
  wv_front <- c(wv_front, max(plants$d))
  
  # If tracking per-seed dispersal contributions, get contributions for foremost individual
  # Then purge ID list of seeds/plants that aren't in the plant or seed bank data frames
  if(sets_dCnts == TRUE){
    wv_dCnts <- rbind(wv_dCnts, wv_ids[wv_ids$ID == plants$ID[which.max(plants$d)], ])
    wv_ids <- wv_ids[wv_ids$ID %in% c(plants$ID, seedsSB$ID), ]
    remove(wv_ids_wd, wv_ids_ad)}
  
  # Plot density over space
  if(i == 1){
    wv_plots <- list()}
  if(sets_plot == TRUE){
    wv_plots[[i]] <- plants$d}
  gc()
  
  # Print procedure progress
  shell("cls")
  cat(paste0(i - 1, "/", sets_years, " (", round(i/(sets_years + 1), 3)*100, "%) complete with ",
             round(as.numeric(difftime(Sys.time(), wv_time, units = "mins")), 2),
             " minutes elapsed"))}

# Close cluster
stopCluster(cl)

# Remove variables no longer in use after run is completed
remove(plants, seedsEst, seedsNew, seedsSB, wv_ids, i, nests, nestsR)

# Print final completion message
shell("cls")
cat(paste0("Procedure complete with ",
           round(as.numeric(difftime(Sys.time(), wv_time, units = "mins")), 2),
           " minutes elapsed"))

