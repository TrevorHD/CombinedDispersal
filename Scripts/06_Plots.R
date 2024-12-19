##### [F2] Create function to plot parameter elasticity ---------------------------------------------------

# Calculate elasticity for each parameter under a 20% decrease
temp1 <- as.numeric(wv_elas[1, ])
temp1 <- ((temp1 - temp1[1])/temp1[1])/(-0.2)
temp1 <- (temp1/sum(abs(temp1)))[-1]

# Re-arrange vector order to match parameter order below
temp1 <- c(temp1[4:6], temp1[1:3], temp1[7:14])
temp1 <- data.frame(var = c("Prob. of seed removal by ants", "Prob. seed predation if removed by ants",
                            "Prob. of seed release",  "Intercept flower head height",
                            "Size-slope flower head height", "SD flower head height",
                            "Surrounding vegetation height", "Mean seed terminal velocity",
                            "SD seed terminal velocity", "Mean wind speed", "SD wind speed",
                            "Maximum thistle density", "Ant nest density", "Maximum ant search radius"),
                    vals = rev(temp1))

# Set plotting colours similar to those in previous figure
temp2 <- c(rep("green", 3), rep("green4", 8), "purple", rep("green4", 2))

# Prepare graphics device
tiff(filename = "Figures/Figure 2.tif", width = 2000, height = 1500, units = "px",
     res = 800, compression = "lzw")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(1500, 2000)
pushViewport(viewport(layout = gly))

# Plot elasticity for each parameter
pushViewport(vp = viewport(layout.pos.row = 1:1500, layout.pos.col = 1:2000))
par(fig = gridFIG())
par(new = TRUE, mar = c(1.5, 4.7, 0.2, 0.5), cex = 0.4, cex.axis = 0.75)
barplot(temp1$vals, space = 0.7, horiz = TRUE, mgp = c(0.3, -0.25, 0), col = rev(temp2),
        xlab = "Scaled Elasticity", tcl = -0.15, xlim = c(-0.4, 0.2))
axis(2, at = seq(1.35, 23.30, length.out = 14), labels = rev(temp1$var),
     las = 2, mgp = c(0.85, 0.2, 0), tcl = 0)
abline(v = 0)
box()
popViewport()

# Deactivate grid layout; finalise graphics save
dev.off()

# Remove plotting objects that are no longer needed
remove(gly, temp1)





##### [F3] Plot wavespeeds over a span of parameter values ------------------------------------------------

# Prepare graphics device
tiff(filename = "Figures/Figure 3.tif", width = 3540, height = 2140, units = "px",
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
plot.span(pert_prob, wv_pert$wv_pert_4, left = TRUE, bottom = FALSE, type = "prob")
popViewport()

# Plot wavespeed when scaling prob. of predation if removed by ants
pushViewport(vp = viewport(layout.pos.row = 741:1440, layout.pos.col = 1:700))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_prob, wv_pert$wv_pert_5, left = TRUE, bottom = FALSE, type = "prob")
popViewport()

# Plot wavespeed when scaling prob. of seed release from capitulum
pushViewport(vp = viewport(layout.pos.row = 1441:2140, layout.pos.col = 1:700))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_prob, wv_pert$wv_pert_6, left = TRUE, bottom = TRUE, type = "prob")
popViewport()

# Plot wavespeed when scaling head height intercept
pushViewport(vp = viewport(layout.pos.row = 40:740, layout.pos.col = 701:1400))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_1, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling head height size-slope
pushViewport(vp = viewport(layout.pos.row = 741:1440, layout.pos.col = 701:1400))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_2, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling head height SD
pushViewport(vp = viewport(layout.pos.row = 1441:2140, layout.pos.col = 701:1400))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_3, left = FALSE, bottom = TRUE, type = "real")
popViewport()

# Plot wavespeed when scaling surrounding vegetation height
pushViewport(vp = viewport(layout.pos.row = 40:740, layout.pos.col = 1401:2100))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_7, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling mean seed terminal velocity
pushViewport(vp = viewport(layout.pos.row = 741:1440, layout.pos.col = 1401:2100))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_8, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling SD seed terminal velocity
pushViewport(vp = viewport(layout.pos.row = 1441:2140, layout.pos.col = 1401:2100))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_9, left = FALSE, bottom = TRUE, type = "real")
popViewport()

# Plot wavespeed when scaling mean wind speed
pushViewport(vp = viewport(layout.pos.row = 40:740, layout.pos.col = 2101:2800))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_10, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling SD wind speed
pushViewport(vp = viewport(layout.pos.row = 741:1440, layout.pos.col = 2101:2800))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_11, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling max thistle density
pushViewport(vp = viewport(layout.pos.row = 1441:2140, layout.pos.col = 2101:2800))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_12, left = FALSE, bottom = TRUE, type = "dens")
popViewport()

# Plot wavespeed when scaling ant nest density 
pushViewport(vp = viewport(layout.pos.row = 40:740, layout.pos.col = 2801:3500))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_13, left = FALSE, bottom = FALSE, type = "real")
popViewport()

# Plot wavespeed when scaling max detection range from ant nests
pushViewport(vp = viewport(layout.pos.row = 741:1440, layout.pos.col = 2801:3500))
par(fig = gridFIG())
par(new = TRUE)
plot.span(pert_real, wv_pert$wv_pert_14, left = FALSE, bottom = TRUE, type = "real")
popViewport()

# Get correlation coefficients and p-values
temp1 <- c()
temp2 <- c()
temp3 <- c(rep("pert_prob", 3), rep("pert_real", 11))
temp4 <- c(4:6, 1:3, 7:14)
for(i in 1:ncol(wv_pert)){
    temp1 <- append(temp1, cor.test(get(temp3[i]), wv_pert[, temp4[i]])$estimate)
    temp2 <- append(temp2, cor.test(get(temp3[i]), wv_pert[, temp4[i]])$p.value)}
temp1 <- paste("r =", sprintf("%.3f", round(temp1, 3)))
for(i in 1:ncol(wv_pert)){
    temp2[i] <- ifelse(as.numeric(temp2[i]) < 0.001, "p < 0.001",
                       paste("p =", sprintf("%.3f", round(as.numeric(temp2[i]), 3))))}

# Add plot labels to reference in figure caption
temp5 <- c(rep(0.189, 3), rep(0.387, 3), rep(0.585, 3), rep(0.783, 3), rep(0.980, 2))
temp6 <- c(rep(c(0.966, 0.638, 0.310), times = 4), 0.966, 0.638)
grid.text(letters[1:14], x = temp5, y = temp6, just = "right", gp = gpar(fontsize = 5))

# Add correlation coefficients and p-values
temp7 <- c(rep(0.069, 3), rep(0.267, 3), rep(0.465, 3), rep(0.663, 3), rep(0.861, 2))
temp8 <- c(rep(c(0.966, 0.638, 0.310), times = 4), 0.966, 0.638)
grid.text(temp1, x = temp7, y = temp8, gp = gpar(fontsize = 4))
temp9 <- c(rep(0.069, 3), rep(0.267, 3), rep(0.465, 3), rep(0.663, 3), rep(0.861, 2))
temp10 <- c(rep(c(0.946, 0.618, 0.290), times = 4), 0.946, 0.618)
grid.text(temp2, x = temp9, y = temp10, gp = gpar(fontsize = 4))

# Deactivate grid layout; finalise graphics save
dev.off()

# Remove plotting objects that are no longer needed
remove(gly, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, i)





##### [FS1, FS2, FS3] Plot snapshots, GIFs of wave movement -----------------------------------------------

# Load simulation settings before running any wavespeed code
sets_plot <- TRUE      # (De)Activate wave plotting
sets_trim <- TRUE      # (De)Activate core area trimming
sets_years <- 100      # Number of years to simulate
sets_trimD <- 1000     # Distance (m) behind wavefront to trim
sets_dCnts <- FALSE    # Track wind/ant dispersal contributions for each seed

# Set demography/dispersal parameters not covered in ADSP or WDSP
misc_tDens <- 10       # Max thistle density per metre
misc_nDens <- 0.2      # Ant nest density (nests/m)
misc_nRange <- 25      # Max detection range (m) from ant nests

# Set seed and demographic/dispersal parameters
set.seed(92633)
aVec <- adsp.param(aNum = 1, aVal = 1)
wVec <- wdsp.param(wNum = 1, wVal = 1)

# Run 04_Wavespeeds before running the code below

# Prepare graphics device
tiff(filename = "Figures/Figure S1.tif", width = 2200, height = 2700, units = "px",
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
plot.wave(type = "line", snap = TRUE, snap_t = 10)
popViewport()

# Plot wavefront for t=20
pushViewport(vp = viewport(layout.pos.row = 500:970, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(0.6, 1.8, 0.0, 0.5))
plot.wave(type = "line", snap = TRUE, snap_t = 20)
popViewport()

# Plot wavefront for t=30
pushViewport(vp = viewport(layout.pos.row = 975:1445, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(0.6, 1.8, 0.0, 0.5))
plot.wave(type = "line", snap = TRUE, snap_t = 30)
popViewport()

# Plot wavefront for t=40
pushViewport(vp = viewport(layout.pos.row = 1450:2000, layout.pos.col = 1:600))
par(fig = gridFIG())
par(new = TRUE, mar = c(1.5, 1.8, 0.0, 0.5))
plot.wave(type = "line", snap = TRUE, snap_t = 40, snap_bottom = TRUE)
popViewport()

# Deactivate grid layout; finalise graphics save
dev.off()

# Remove plotting objects that are no longer needed
remove(gly)

# Generate GIFs of population spread (will be included in repo, but not publication)
save_gif(plot.wave("hist"), "Figures/Figure S2.gif", delay = 0.2, width = 1280, height = 720, res = 144)
save_gif(plot.wave("line"), "Figures/Figure S3.gif", delay = 0.2, width = 1280, height = 720, res = 144)

