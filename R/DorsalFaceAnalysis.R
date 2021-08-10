### Dorsal Face Analysis ###


### Perform GPA alignment for all three datasets ###
DorsalFace_20_GPA <- gpagen(Croc_DorsalFace_SubCurve_20,
                            curves=Curves_20pt_sliders,
                            ProcD=TRUE)
###

### Perform PCA for all three datasets ###
DorsalFace_20_PCA <- gm.prcomp(DorsalFace_20_GPA$coords) #snout profile only dataset PCA
###

### Test correclation between PCs and centroid size ###
cor(DorsalFace_20_PCA$x,log(Croc_Lateral_GPA$Csize))

p.cor <- NULL
for (i in 1:ncol(DorsalFace_20_PCA$x)){
  p.cor <- rbind(p.cor,cor.test(log(Croc_Lateral_GPA$Csize),DorsalFace_20_PCA$x[,i])$p.value)
}
###

### Test for mean shape differences among ecomorphs at each ontogenetic stage ###
DorsalFace_20_Mean_shape_LMs <- list()
# for loop to iterate through each ontogenetic class and perform P-ANOVA for ecomorph differences in mean shape
for (i in 1:length(Ontogeny_list)){ 
  temp_stage <- Ontogeny_list[i]
  temp_grep <- grep(temp_stage, lateral_covariates$Ontogeny)
  temp_lm <- procD.lm(DorsalFace_20_GPA$coords[,,temp_grep] ~ lateral_covariates$simple_ecomorph[temp_grep])
  
  DorsalFace_20_Mean_shape_LMs[[temp_stage]] <- temp_lm
}
#
lapply(DorsalFace_20_Mean_shape_LMs, "summary") #Summary of differences in mean shape among ecomorphs at each individual stage class
p.adjust(c(0.065,0.001,0.001,0.001,0.001,0.002)) #Correction for multiple comparisons
###

### Test for differences in ontogenetic trajectories of ecomorphs & species ###
# Whole ontogeny, across ecomorphs
DorsalFace_Ecomorph_Allometry_LM <- procD.lm(DorsalFace_20_GPA$coords ~ log(Croc_Lateral_GPA$Csize) * lateral_covariates$simple_ecomorph)
summary(DorsalFace_Ecomorph_Allometry_LM)
summary(pairwise(DorsalFace_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph, covariate = log(Croc_Lateral_GPA$Csize)), test.type = "VC")
DorsalFace_20_CACs <- plotAllometry(DorsalFace_Ecomorph_Allometry_LM, log(Croc_Lateral_GPA$Csize), method="CAC")$CAC #calculation of CAC
#
# Whole ontogeny, across species
DorsalFace_Species_Allometry_LM <- procD.lm(DorsalFace_20_GPA$coords ~ log(Croc_Lateral_GPA$Csize) * lateral_covariates$species)
summary(DorsalFace_Species_Allometry_LM)
summary(pairwise(DorsalFace_Species_Allometry_LM, groups = lateral_covariates$species, covariate = log(Croc_Lateral_GPA$Csize)), test.type = "VC")
#
# Embryonic development only, across ecomorphs
Embryos_DorsalFace_Ecomorph_Allometry_LM <- procD.lm(DorsalFace_20_GPA$coords[,,Embryos_list] ~ log(Croc_Lateral_GPA$Csize[Embryos_list]) * lateral_covariates$simple_ecomorph[Embryos_list])
summary(Embryos_DorsalFace_Ecomorph_Allometry_LM) #P-ANOVA table
summary(pairwise(Embryos_DorsalFace_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[Embryos_list])) #Pairwise differences in mean shape
summary(pairwise(Embryos_DorsalFace_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[Embryos_list], covariate = log(Croc_Lateral_GPA$Csize[Embryos_list])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#
# Embryonic development only, across species
Embryos_DorsalFace_Species_Allometry_LM <- procD.lm(DorsalFace_20_GPA$coords[,,Embryos_list[-29]] ~ log(Croc_Lateral_GPA$Csize[Embryos_list[-29]]) * lateral_covariates$species[Embryos_list[-29]])
summary(Embryos_DorsalFace_Species_Allometry_LM) #P-ANOVA table
# must exclude the Gavilis embryos from pairwise comparisons, as with only a single individual the comparisons are impossible otherwise
summary(pairwise(Embryos_DorsalFace_Species_Allometry_LM, groups = lateral_covariates$species[Embryos_list[-29]])) #Pairwise differences in mean shape
summary(pairwise(Embryos_DorsalFace_Species_Allometry_LM, groups = lateral_covariates$species[Embryos_list[-29]], covariate = log(Croc_Lateral_GPA$Csize[Embryos_list[-29]])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#
# Post-hatching growth only, across ecomorphs
Posthatching_DorsalFace_Ecomorph_Allometry_LM <- procD.lm(DorsalFace_20_GPA$coords[,,-Embryos_list] ~ log(Croc_Lateral_GPA$Csize[-Embryos_list]) * lateral_covariates$simple_ecomorph[-Embryos_list])
summary(Posthatching_DorsalFace_Ecomorph_Allometry_LM) #P-ANOVA table
summary(pairwise(Posthatching_DorsalFace_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[-Embryos_list])) #Pairwise differences in mean shape
summary(pairwise(Posthatching_DorsalFace_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[-Embryos_list], covariate = log(Croc_Lateral_GPA$Csize[-Embryos_list])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#
# Post-hatching growth only, across species
Posthatching_DorsalFace_Species_Allometry_LM <- procD.lm(DorsalFace_20_GPA$coords[,,-Embryos_list] ~ log(Croc_Lateral_GPA$Csize[-Embryos_list]) * lateral_covariates$species[-Embryos_list])
summary(Posthatching_DorsalFace_Species_Allometry_LM) #P-ANOVA table
summary(pairwise(Posthatching_DorsalFace_Species_Allometry_LM, groups = lateral_covariates$species[-Embryos_list])) #Pairwise differences in mean shape
summary(pairwise(Posthatching_DorsalFace_Species_Allometry_LM, groups = lateral_covariates$species[-Embryos_list], covariate = log(Croc_Lateral_GPA$Csize[-Embryos_list])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#

#correction of p values for pairwise comparisons
p.adjust(c(0.032,0.009,0.004,0.346,0.242,0.083,0.140,0.049,0.027))
p.adjust(c(0.037,0.009,0.019,0.493,0.416,0.128,0.082,0.039,0.040)) #Correction for multiple comparisons among ecomorphs
###


###Calculate slopes for each trajectory#
DorsalFace_Slopes <- list()
# for loop to iterate through each ontogenetic class and perform P-ANOVA for ecomorph differences in mean shape
for(j in 1:3){
  temp_group <- c("Complete", "Embryonic", "Post-hatching")[j]
  temp_CAC <- list(DorsalFace_20_CACs, DorsalFace_20_CACs[Embryos_list], DorsalFace_20_CACs[-Embryos_list])[[j]]
  temp_CS <- list(Croc_Lateral_GPA$Csize, Croc_Lateral_GPA$Csize[Embryos_list], Croc_Lateral_GPA$Csize[-Embryos_list])[[j]]
  temp_cov <- list(lateral_covariates, lateral_covariates[Embryos_list,], lateral_covariates[-Embryos_list,])[[j]]
  slope_matrix <- matrix(data=NA,nrow=3,ncol=2,dimnames = list(c(unique(lateral_covariates$simple_ecomorph)),c("Intercept","Slope")))
  for (i in 1:length(unique(lateral_covariates$simple_ecomorph))){ 
    temp_eco <- unique(temp_cov$simple_ecomorph)[i]
    temp_grep <- grep(temp_eco, temp_cov$simple_ecomorph)
    temp_lm <- lm(temp_CAC[temp_grep] ~ log(temp_CS[temp_grep]))
    slope_matrix[temp_eco,] <- temp_lm$coefficients
  }
  DorsalFace_Slopes[[temp_group]] <- slope_matrix
}
DorsalFace_Slopes
###