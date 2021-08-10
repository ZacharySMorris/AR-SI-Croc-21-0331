### Dependancies ###
library(geomorph)
library(scales)
library(mcp)
###

### Load landmark and covariate data ###
load("data/covariates.RData") #load in the specimen co-variate data
load("data/landmarks.RData") #load in complete skull landmarks
load("data/snoutprofile_curve_landmarks.RData") #load in dorsal face / snout profile curve semi-landmarks
load("data/skulltable_curve_landmarks.RData") #load in skull table curve semi-landmarks
load("data/completelateralskullsliders.RData") #load in complete skull semi-landmark sliding data
load("data/curve20sliders.RData") #load in 20 point curve only semi-landmark sliding data
###

### Create important lists for dividing up ontogeny
# list of all ontogenetic classes
Ontogeny_list <- unique(lateral_covariates$Ontogeny)
#
# list of all embryos
toMatch <- c("Mid-Stage Embryo","Late-Stage Embryo")
Embryos_list <- grep(paste(toMatch,collapse = "|"),lateral_covariates$Ontogeny)
#
###
