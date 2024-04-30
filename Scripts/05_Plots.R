##### [FS1] Plot snapshots of wave movement ----------------------------------------------------------------

# Prepare graphics device
tiff(filename = "Figures/Figure S1.tif", width = 2000, height = 3000, units = "px",
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





##### [FS2] Plot wavespeeds over a span of parameter values ------------------------------------------------

# Set placeholder data
test1 <- c(1:11)

# Prepare graphics device
tiff(filename = "Figures/Figure S2.tif", width = 3540, height = 2140, units = "px",
     res = 800, compression = "lzw")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(2140, 3540)
pushViewport(viewport(layout = gly))

# Plot wavespeed when scaling prob. of seed removal by ants
pushViewport(vp = viewport(layout.pos.row = 40:740, layout.pos.col = 1:700))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multPr, test1, left = TRUE, bottom = FALSE, type = "prob")
popViewport()

# Plot wavespeed when scaling prob. of surviving predation if removed by ants
pushViewport(vp = viewport(layout.pos.row = 741:1440, layout.pos.col = 1:700))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multPr, test1, left = TRUE, bottom = FALSE, type = "prob")
popViewport()

# Plot wavespeed when scaling prob. of seed release from capitulum
pushViewport(vp = viewport(layout.pos.row = 1441:2140, layout.pos.col = 1:700))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multPr, test1, left = TRUE, bottom = TRUE, type = "prob")
popViewport()

# Plot wavespeed when scaling head height intercept
pushViewport(vp = viewport(layout.pos.row = 40:740, layout.pos.col = 701:1400))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling head height size-slope
pushViewport(vp = viewport(layout.pos.row = 741:1440, layout.pos.col = 701:1400))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling head height SD
pushViewport(vp = viewport(layout.pos.row = 1441:2140, layout.pos.col = 701:1400))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = TRUE, type = "real")
popViewport()

# Plot wavespeed when scaling surrounding vegetation height
pushViewport(vp = viewport(layout.pos.row = 40:740, layout.pos.col = 1401:2100))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling mean seed terminal velocity
pushViewport(vp = viewport(layout.pos.row = 741:1440, layout.pos.col = 1401:2100))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling SD seed terminal velocity
pushViewport(vp = viewport(layout.pos.row = 1441:2140, layout.pos.col = 1401:2100))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = TRUE, type = "real")
popViewport()

# Plot wavespeed when scaling mean wind speed
pushViewport(vp = viewport(layout.pos.row = 40:740, layout.pos.col = 2101:2800))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling SD wind speed
pushViewport(vp = viewport(layout.pos.row = 741:1440, layout.pos.col = 2101:2800))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling max thistle density
pushViewport(vp = viewport(layout.pos.row = 1441:2140, layout.pos.col = 2101:2800))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = TRUE, type = "real")
popViewport()

# Plot wavespeed when scaling ant nest density 
pushViewport(vp = viewport(layout.pos.row = 40:740, layout.pos.col = 2801:3500))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling max detection range from ant nests
pushViewport(vp = viewport(layout.pos.row = 741:1440, layout.pos.col = 2801:3500))
par(fig = gridFIG())
par(new = TRUE)
plot.span(multRe, test1, left = FALSE, bottom = TRUE, type = "real")
popViewport()

# Deactivate grid layout; finalise graphics save
dev.off()





##### [FS3] Create function to plot parameter elasticity ---------------------------------------------------

# Create placeholder data with elasticity for each scalable parameter
temp <- data.frame(var = c("Mean initial rosette size", "SD initial rosette size", "Prob. rosette survival",
                           "Intercept rosette size after growth", "Size-slope rosette size after growth",
                           "SD rosette size after growth", "Prob. rosette bolting/flowering",
                           "Intercept number of flower heads", "Size-slope number of flower heads",
                           "SD number of flower heads", "Intercept flower head height",
                           "Size-slope flower head height", "SD flower head height", "Max thistle density",
                           "Seeds per flower head", "Prob. establishment from seed", "Prob. seed survival in seed bank",
                           "Prob. seed establishing from seed bank", "Prob. seed entering seed bank",
                           "Prop. seeds surviving florivory", "Prop. seeds released",
                           "Mean seed terminal velocity", "SD seed terminal velocity", "Mean wind speed",
                           "SD wind speed", "Vegetation height", "Prob. seed removal by ants",
                           "Prob. surviving post-dispersal predation", "Ant nest density",
                           "Maximum ant search radius"),
                   vals = seq(-1, 1, length.out = 30))

# Prepare graphics device
tiff(filename = "Figures/Figure S3.tif", width = 2000, height = 2000, units = "px",
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
par(new = TRUE, mar = c(1.5, 4.7, 0.2, 0.5), cex = 0.4, cex.axis = 0.75)
barplot(temp$vals, space = 0.7, horiz = TRUE, mgp = c(0.3, -0.25, 0),
        xlab = "Elasticity", tcl = -0.15)
axis(2, at = seq(1.35, 50.75, length.out = 30), labels = temp$var,
     las = 2, mgp = c(0.85, 0.2, 0), tcl = 0)
box()
popViewport()

# Deactivate grid layout; finalise graphics save
dev.off()





##### [FS4, FS5] Plot moving wave over time ---------------------------------------------------------------

# Generate GIFs of population spread (will be included in repo, but not publication)
save_gif(plot.wave("hist"), "Figures/Figure S4.gif", delay = 0.2, width = 1280, height = 720, res = 144)
save_gif(plot.wave("density"), "Figures/Figure S5.gif", delay = 0.2, width = 1280, height = 720, res = 144)

