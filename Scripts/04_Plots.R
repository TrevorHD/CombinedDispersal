##### Create function to plot wavefront over time ---------------------------------------------------------

# Function to generate GIF of population spread (nYear = 100 recommended)
generatePlots <- function(type = "hist", snap = FALSE, snap_t = NULL, snap_bottom = FALSE){
  
  # Set x-axis maximum based on wavefront at end of simulation (or a specified number if using snapshots)
  # Set other graphical parameters based on whether or not there will be multiple adjacent plots
  if(snap == FALSE){
    max_scale <- (ceiling(max(wv_plots[[length(wv_plots)]])/1000) + 1)*1000}
  if(snap == TRUE){
    max_scale <- 3000}
  xlab <- ifelse(snap_bottom == TRUE | snap == FALSE, "Distance (m)", "")
  tsize <- ifelse(snap == FALSE, 1, 0.5)
  
  if(snap == FALSE){
    plotSeq <- 1:length(wv_plots)}
  if(snap == TRUE){
    plotSeq <- snap_t}
  
  # Generate plot using position data at each timestep
  for(i in plotSeq){
    
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
    
    # Plot wave movement as either histogram or density relative to maximum
    if(type == "hist"){
      hist(c(lower, positions), breaks = seq(0, max_scale, by = 20), xlim = c(0, max_scale), ylim = c(0, 215),
           xaxt = "n", yaxt = "n", xlab = NULL, ylab = "Count", main = "")
      axis(2, at = seq(0, 200, 50))}
    if(type == "density"){
      plotdata <- hist(c(lower, positions), breaks = seq(0, max_scale, by = 20), plot = FALSE)
      plot(plotdata$mids, plotdata$density/max(plotdata$density),
           type = "l", xlim = c(0, max_scale), ylim = c(0, 1.075),
           xaxt = "n", yaxt = "n", xlab = "", ylab = "Relative Density", main = "")
      axis(2, at = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"))}
    if(snap == FALSE){
      axis(1, at = seq(0, max_scale, 1000))
      title(xlab = xlab)}
    if(snap == TRUE){
      if(snap_bottom == FALSE){
        axis(1, at = seq(0, max_scale, 1000), labels = FALSE)}
      if(snap_bottom == TRUE){
        axis(1, at = seq(0, max_scale, 1000), mgp = c(0, -0.1, 0))
        title(xlab = xlab, mgp = c(0.5, 0, 0))}}
    if(type == "hist"){
      text(x = max_scale*0.95, y = 200, paste0("t = ", i), cex = tsize)
      box()}
    if(type == "density"){
      text(x = max_scale*0.95, y = 1, paste0("t = ", i), cex = tsize)}}}

    



##### [F1] Plot snapshots of wave movement ---------------------------------------------------------

# Prepare graphics device
tiff(filename = "Figure 1.tif", width = 2000, height = 3000, units = "px",
     res = 800, compression = "lzw")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(2000, 600)
pushViewport(viewport(layout = gly))

# CN non-warmed: mean vs distribution
pushViewport(vp = viewport(layout.pos.row = 1:500, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(0.6, 1.8, 0.4, 0.5), cex.axis = 0.45, cex.lab = 0.5,
    tcl = -0.15, mgp = c(0.85, 0.2, 0))
generatePlots(type = "density", snap = TRUE, snap_t = 10)
popViewport()

# CN warmed: mean vs distribution
pushViewport(vp = viewport(layout.pos.row = 500:970, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(0.6, 1.8, 0.0, 0.5))
generatePlots(type = "density", snap = TRUE, snap_t = 20)
popViewport()

# CA non-warmed: mean vs distribution
pushViewport(vp = viewport(layout.pos.row = 975:1445, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(0.6, 1.8, 0.0, 0.5))
generatePlots(type = "density", snap = TRUE, snap_t = 30)
popViewport()

# CA warmed: mean vs distribution
pushViewport(vp = viewport(layout.pos.row = 1450:2000, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(1.5, 1.8, 0.0, 0.5))
generatePlots(type = "density", snap = TRUE, snap_t = 40, snap_bottom = TRUE)
popViewport()

# Deactivate grid layout; finalise graphics save
dev.off()

# Generate GIFs of population spread (will be included in repo, but not publication)
save_gif(generatePlots("hist"), "Spread1.gif", delay = 0.2, width = 1280, height = 720, res = 144)
save_gif(generatePlots("density"), "Spread2.gif", delay = 0.2, width = 1280, height = 720, res = 144)





##### Create function to plot parameter elasticity --------------------------------------------------------

# Create placeholder data with elasticity for each scalable parameter
temp <- data.frame(var = c("Vegetation height", "Mean wind speed", "SD wind speed",
                           "Mean terminal velocity", "SD terminal velocity",
                           "Intercept prob. ant dispersal", "Slope prob. ant dispersal",
                           "Seeds per flower head", "Prob. estab. from seed",
                           "Prob. seed predation", "Prob. seed entering SB",
                           "Prob. estab. from SB", "Prob. survival in SB",
                           "Mean rosette size", "SD rosette size", "Prob. flowering",
                           "Intercept num. heads", "Slope num. heads", "Mean head height",
                           "SD head height", "Prob. survival"),
                   vals = seq(-1, 1, length.out = 21))

# Prepare graphics device
tiff(filename = "Figure 2.tif", width = 2000, height = 2000, units = "px",
     res = 800, compression = "lzw")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(2000, 2000)
pushViewport(viewport(layout = gly))

# Plot elasticity values
pushViewport(vp = viewport(layout.pos.row = 1:2000, layout.pos.col = 1:2000))
par(fig = gridFIG())
par(new = TRUE, mar = c(1.5, 3.7, 0.2, 0.5), cex = 0.4, cex.axis = 0.75)
barplot(temp$vals, space = 0.7, horiz = TRUE, mgp = c(0.3, -0.25, 0),
        xlab = "Elasticity", tcl = -0.15)
axis(2, at = seq(1.25, 35.3, length.out = 21), labels = temp$var,
     las = 2, mgp = c(0.85, 0.2, 0), tcl = 0)
box()
popViewport()

# Deactivate grid layout; finalise graphics save
dev.off()

