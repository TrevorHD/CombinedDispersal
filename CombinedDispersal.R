# Randomly generate ant nests on the landscape
# Place single mature individual at centre
# Individual reproduces and then seeds disperse
# For all dispersed seeds, model dispersal with a predation coefficient
# Run demo matrix to get germination rate and growth
# Repeat previous 3 steps


library(SuppDists)
library(MASS)
library(tidyverse)
library(truncnorm)

#plants <- data.frame(ID = 1, r = 0, theta = 0, x = 0, y = 0, size = 0.8, stage = 2)
plants <- data.frame(r = c(1, 4, 6, 2), theta = c(pi/2, pi, 3*pi/2, 2*pi),
                     x = c(0, -4, 0, 2), y = c(1, 0, -6, 0), stage = c(1, 2, 2, 2))

seeds <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(seeds) <- c("r", "theta", "x", "y", "germ")

# Placeholder dispersal kernel
kern <- function(n, x = 0, y = 0){
  theta <- sample(seq(0.1, 2*pi, by = pi/100), n, replace = TRUE)
  r <- rlnorm(n)
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







# 11 points
# 10 bins
seq(0, 1.5, by = 0.15)
# size is integer between 1 and 10

# two classes: rosette (y1) and adult (y2)


# Demography
# In a given iteration, run reproduction, then germination, then survival, then growth
demo <- function(dType){
  if(dType == "reproduction"){
    return(100)}
  if(dType == "germination"){
    prob <- 0.05
    outcome <- sample(c(0, 1), 1, prob = c(1 - prob, prob))
    return(outcome)}
  if(dType == "survival"){
    prob <- 0.92
    outcome <- sample(c(0, 1), 1, prob = c(1 - prob, prob))
    return(outcome)}
  if(dType == "growth"){
    return(rtruncnorm(1, a = 0.4, b = 1.7, mean = 1))}}





