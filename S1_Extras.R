##### 2D expansion ----------------------------------------------------------------------------------------

# Estimate dispersal distances from given point; assume 1m plant height
kern <- function(n, x0 = 0, y0 = 0){
  r <- WALD.b(n, 1)
  theta <- sample(seq(0.1, 2*pi, by = pi/100), n, replace = TRUE)
  x <- r*cos(theta) + x0
  y <- r*sin(theta) + y0
  return(data.frame(theta, r, x, y))}

# Initialise data; start with a single rosette and adult
plants <- data.frame(r = c(0.01, 0.01), theta = c(0, 0), x = c(0.01, 0.01), y = c(0.01, 0.01), stage = c(0, 2))
seeds <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(seeds) <- c("r", "theta", "x", "y", "germ")
vals <- c()

# Placeholder nest generation
# Max range = 100 m
#nests <- function(n){
#  theta <- sample(seq(0.1, 2*pi, by = pi/100), n, replace = TRUE)
#  r <- sample(seq(0.1, 100, by = 0.1), n, replace = TRUE)
#  x <- r*cos(theta)
#  y <- r*sin(theta)
#  return(data.frame(theta, r, x, y))}

# Simulate dispersal and seed survival
for(i in 1:nrow(plants)){
  if(plants$stage[i] == 2){
    n <- demo("reproduction")
    newSeeds <- kern(n, plants$x[i], plants$y[i])
    newSeeds <- cbind(newSeeds, rep(1, n))
    names(newSeeds)[5] <- "germ"
    newSeeds <- newSeeds[newSeeds$germ == 1, ]
    seeds <- rbind(seeds, newSeeds)}}

# Kill adults after they reproduce
plants <- plants[plants$stage != 2, ]

# Rosettes from previous year become adults
plants$stage[plants$stage == 0] <- 2

# Surviving seeds become rosettes
seeds$stage <- rep(0, nrow(seeds))
seeds <- seeds[, !names(seeds) == c("germ")]
plants <- rbind(plants, seeds)

# Kill rosettes if density is too high
# Do this by sorting so that adults come first and are prioritised
plants %>% 
  mutate(xbin = as.integer(cut(x, breaks = seq(floor(min(plants$x)), ceiling(max(plants$x)), 1))),
         ybin = as.integer(cut(y, breaks = seq(floor(min(plants$y)), ceiling(max(plants$y)), 1)))) %>% 
  arrange(xbin, ybin, desc(stage)) %>% 
  group_by(xbin, ybin) %>%
  slice(1:pmin(10, n())) %>% 
  data.frame() -> plants
plants <- plants[, !names(plants) %in% c("xbin", "ybin")]

# Reset seeds
# Can change this later to account for seed bank
seeds <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(seeds) <- c("r", "theta", "x", "y", "germ")

# Plot density over space
# Note: r is relative to last point of dispersal, so we calculate distance manually
plot(plants$x, plants$y, xlim = c(-200, 200), ylim = c(-200, 200),
     col = rgb(r = 0, g = 0, b = 0, alpha = 0.2), pch = 16)
vals <- c(vals, max(sqrt((plants$x)^2 + (plants$y)^2)))

# Get wavespeeds
diff(vals)

