##### [F1] Plot snapshots of wave movement ----------------------------------------------------------------

# Prepare graphics device
tiff(filename = "Figures/Figure 1.tif", width = 2000, height = 3000, units = "px",
     res = 800, compression = "lzw")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(2000, 600)
pushViewport(viewport(layout = gly))

# Plot wavefront for t=10
pushViewport(vp = viewport(layout.pos.row = 1:500, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(0.6, 1.8, 0.4, 0.5), cex.axis = 0.45, cex.lab = 0.5,
    tcl = -0.15, mgp = c(0.85, 0.2, 0))
generatePlots(type = "density", snap = TRUE, snap_t = 10)
popViewport()

# Plot wavefront for t=20
pushViewport(vp = viewport(layout.pos.row = 500:970, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(0.6, 1.8, 0.0, 0.5))
generatePlots(type = "density", snap = TRUE, snap_t = 20)
popViewport()

# Plot wavefront for t=30
pushViewport(vp = viewport(layout.pos.row = 975:1445, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(0.6, 1.8, 0.0, 0.5))
generatePlots(type = "density", snap = TRUE, snap_t = 30)
popViewport()

# Plot wavefront for t=40
pushViewport(vp = viewport(layout.pos.row = 1450:2000, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(1.5, 1.8, 0.0, 0.5))
generatePlots(type = "density", snap = TRUE, snap_t = 40, snap_bottom = TRUE)
popViewport()

# Deactivate grid layout; finalise graphics save
dev.off()





##### [F2] Create function to plot parameter elasticity ---------------------------------------------------

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
tiff(filename = "Figures/Figure 2.tif", width = 2000, height = 2000, units = "px",
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





##### [FS1, FS2] Plot moving wave over time ---------------------------------------------------------------

# Generate GIFs of population spread (will be included in repo, but not publication)
save_gif(generatePlots("hist"), "Figures/Figure S1.gif", delay = 0.2, width = 1280, height = 720, res = 144)
save_gif(generatePlots("density"), "Figures/Figure S2.gif", delay = 0.2, width = 1280, height = 720, res = 144)

