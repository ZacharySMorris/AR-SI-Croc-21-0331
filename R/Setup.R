### Load landmark and covariate data ###
load("data/covariates.RData") #load in the specimen co-variate data
load("data/landmarks.RData") #load in complete skull landmarks
load("data/snoutprofile_curve_landmarks.RData") #load in dorsal face / snout profile curve semi-landmarks
load("data/skulltable_curve_landmarks.RData") #load in skull table curve semi-landmarks
load("data/completelateralskullsliders.RData") #load in complete skull semi-landmark sliding data
load("data/curve20sliders.RData") #load in 20 point curve only semi-landmark sliding data
###

### Perform GPA alignment for all three datasets ###
Croc_Lateral_GPA <- gpagen(estimate.missing(Compiled_Lateral_LM, method = "TPS"),
                           curves=complete.lateral.skull.sliders,
                           ProcD=FALSE)

DorsalFace_20_GPA <- gpagen(Croc_DorsalFace_SubCurve_20,
                            curves=Curves_20pt_sliders,
                            ProcD=TRUE)


SkullTable_20_GPA <- gpagen(Croc_SkullTable_SubCurve_20,
                            curves=Curves_20pt_sliders,
                            ProcD=TRUE)
###

### Perform PCA for all three datasets ###
Croc_Lateral_PCA <- gm.prcomp(Croc_Lateral_GPA$coords) #complete skull dataset PCA

DorsalFace_20_PCA <- gm.prcomp(DorsalFace_20_GPA$coords) #snout profile only dataset PCA

SkullTable_20_PCA <- gm.prcomp(SkullTable_20_GPA$coords) #skull table only dataset PCA
###