##### Create animation of wave movement -------------------------------------------------------------------

# Function to generate GIF of population spread (nYear = 100 recommended)
generatePlots <- function(type = "hist"){
  
  # Set x-axis maximum based on wavefront at end of simulation
  max_scale <- (ceiling(max(wv_plots[[length(wv_plots)]])/1000) + 1)*1000
  
  # Generate plot using position data at each timestep
  for(i in 1:length(wv_plots)){
    
    # Initialise vectors
    lower <- c()
    positions <- wv_plots[[i]]
    if(length(positions) == 0){
      positions <- 0}
    
    # Assume max density behind trimmed wavefront
    # Trim condition here should match that used in the simulation
    if(max(positions > 1000)){
      minBin <- sum(positions > floor(min(positions)) & positions < ceiling(min(positions)))
      lower <- c(rep(0:(floor(min(positions) - 1)), 10) + 0.01,
                 rep(floor(min(positions)) + 0.01, 10 - minBin))
      if(min(positions) < 1){
        lower <- sort(lower)[-c(1:20)]}}
    
    # Plot wave movement as histogram
    if(type == "hist"){
      hist(c(lower, positions), breaks = seq(0, max_scale, by = 20), xlim = c(0, max_scale), ylim = c(0, 250),
           xaxt = "n", yaxt = "n", xlab = "Distance (m)", ylab = "Count", main = "")
      axis(1, at = seq(0, max_scale, 1000))
      axis(2, at = seq(0, 500, 50))
      text(x = max_scale*0.95, y = 230, paste0("t = ", i))
      box()}
    
    # Plot wave movement as density relative to maximum
    if(type == "density"){
      plotdata <- hist(c(lower, positions), breaks = seq(0, max_scale, by = 20), plot = FALSE)
      plot(plotdata$mids, plotdata$density/max(plotdata$density),
           type = "l", xlim = c(0, max_scale), ylim = c(0, 1.25),
           xaxt = "n", yaxt = "n", xlab = "Distance (m)", ylab = "Relative Density", main = "")
      axis(1, at = seq(0, max_scale, 1000))
      axis(2, at = seq(0, 1.25, 0.25))
      text(x = max_scale*0.95, y = 1.15, paste0("t = ", i))}}}

# Generate GIFs of population spread
save_gif(generatePlots("hist"), "Spread1.gif", delay = 0.2, width = 1280, height = 720, res = 144)
save_gif(generatePlots("density"), "Spread2.gif", delay = 0.2, width = 1280, height = 720, res = 144)

