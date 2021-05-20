##### Load libraries and initialise data ------------------------------------------------------------------

# Load libraries
library(SuppDists)
library(MASS)
library(tidyverse)
library(truncnorm)





##### 2D expansion ----------------------------------------------------------------------------------------

plants <- data.frame(r = 0, theta = 0, x = 0, y = 0, stage = 2)

seeds <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(seeds) <- c("r", "theta", "x", "y", "germ")

# Placeholder dispersal kernel
kern <- function(n, x = 0, y = 0){
  theta <- sample(seq(0.1, 2*pi, by = pi/100), n, replace = TRUE)
  r <- rlnorm(n, meanlog = 0.25, sdlog = 1.1)
  x <- r*cos(theta) + x
  y <- r*sin(theta) + y
  return(data.frame(theta, r, x, y))}

# Placeholder nest generation
# Max range = 100 m
nests <- function(n){
  theta <- sample(seq(0.1, 2*pi, by = pi/100), n, replace = TRUE)
  r <- sample(seq(0.1, 100, by = 0.1), n, replace = TRUE)
  x <- r*cos(theta)
  y <- r*sin(theta)
  return(data.frame(theta, r, x, y))}

# for below, separate rosettes and adults first, then rejoin after performing procedures on adults

# Simulate dispersal and seed survival
for(i in 1:nrow(plants)){
  if(plants$stage[i] == 2){
    n <- demo("reproduction")
    newSeeds <- kern(n, plants$x[i], plants$y[i])
    newSeeds <- cbind(newSeeds, replicate(n, demo("germination")))
    names(newSeeds)[5] <- "germ"
    newSeeds <- newSeeds[newSeeds$germ == 1, ]
    seeds <- rbind(seeds, newSeeds)}}

# for below, separate rosettes and adults first, then rejoin after performing procedures on adults

# Kill rosettes if density is too high
for(i in 1:nrow(plants)){
  if(plants$stage[i] == 1){
    dists <- (plants$x[i] - plants$x)^2 + (plants$y[i] - plants$y)^2
    if(length(dists[dists < 1]) > 8){
      plants <- plants[-i, ]
    }
  }
}

# Adults die after reproducing
plants <- plants[plants$stage != 2, ]

# Rosettes from previous year become adults
plants$stage[plants$stage == 1] <- 2

# Surviving seeds become rosettes
seeds <- seeds[, !names(seeds) == c("germ")]
seeds$stage <- rep(1, nrow(seeds))
plants <- rbind(plants, seeds)

# Reset seeds
# Can change this later to account for seed bank
seeds <- data.frame(matrix(ncol = 5, nrow = 0))

plot(plants$x, plants$y, xlim = c(-50, 50), ylim = c(-50, 50), col = rgb(r = 0, g = 0, b = 0, alpha = 0.2), pch = 16)





##### 1D expansion ----------------------------------------------------------------------------------------

# Placeholder dispersal kernel until WALD is implemented
kern <- function(n, d0 = 0){
  d <- rlnorm(n, meanlog = 0.25, sdlog = 1.1) + d0
  return(d)}

# Initialise data; start with a single rosette and adult
plants <- data.frame(d = c(0.01, 0.01), stage = c(1, 2))
seeds <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(seeds) <- c("d", "germ")
vals <- c()

# Simulate dispersal and seed survival
for(i in 1:nrow(plants)){
  if(plants$stage[i] == 2){
    n <- demo("reproduction")
    newSeeds <- data.frame(kern(n, plants$d[i]))
    newSeeds <- cbind(newSeeds, rep(1, n))
    names(newSeeds) <- c("d", "germ")
    newSeeds <- newSeeds[newSeeds$germ == 1, ]
    seeds <- rbind(seeds, newSeeds)}}

# Kill adults after they reproduce
plants <- plants[plants$stage != 2, ]

# Rosettes from previous year become adults
plants$stage[plants$stage == 1] <- 2

# Surviving seeds become rosettes
seeds$stage <- rep(1, nrow(seeds))
seeds <- seeds[, !names(seeds) == c("germ")]
plants <- rbind(plants, seeds)

# Kill rosettes if density is too high
# Do this by sorting so that adults come first and are prioritised
plants %>% 
  mutate(bin = as.integer(cut(d, breaks = seq(0, ceiling(max(plants$d)), 1)))) %>% 
  arrange(bin, desc(stage)) %>% 
  group_by(bin) %>%
  slice(1:pmin(10, n())) %>% 
  data.frame() -> plants
plants <- plants[, !names(plants) == c("bin")]

# Reset seeds
# Can change this later to account for seed bank
seeds <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(seeds) <- c("d", "germ")

# Plot density over space
hist(plants$d, breaks = seq(0, 1000, by = 1), ylim = c(0, 12), xlab = "Distance", ylab = "Density")
vals <- c(vals, max(plants$d))

# Get wavespeeds
diff(vals)











# Demography
# In a given iteration, run reproduction, then germination, then survival, then growth
demo <- function(dType){
  if(dType == "reproduction"){
    nSeed <- 1000
    pPred <- 0.95
    return(nSeed - nSeed*pPred)}
  if(dType == "survival"){
    prob <- 0.95
    outcome <- sample(c(0, 1), 1, prob = c(1 - prob, prob))
    return(outcome)}
  if(dType == "growth"){
    return(rtruncnorm(1, a = 0.4, b = 1.7, mean = 1))}}





