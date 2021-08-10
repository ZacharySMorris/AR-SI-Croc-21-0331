### Complete Skull Analysis ###


### Perform GPA alignment for all three datasets ###
Croc_Lateral_GPA <- gpagen(estimate.missing(Compiled_Lateral_LM, method = "TPS"),
                           curves=complete.lateral.skull.sliders,
                           ProcD=FALSE)
###

### Perform PCA for all three datasets ###
Croc_Lateral_PCA <- gm.prcomp(Croc_Lateral_GPA$coords) #complete skull dataset PCA
###

### Test correclation between PCs and centroid size ###
cor(Croc_Lateral_PCA$x,log(Croc_Lateral_GPA$Csize))

p.cor <- NULL
for (i in 1:ncol(Croc_Lateral_PCA$x)){
  p.cor <- rbind(p.cor,cor.test(log(Croc_Lateral_GPA$Csize),Croc_Lateral_PCA$x[,i])$p.value)
}
###

### Test for mean shape differences among ecomorphs at each ontogenetic stage ###
Croc_Lateral_Mean_shape_LMs <- list()
# for loop to iterate through each ontogenetic class and perform P-ANOVA for ecomorph differences in mean shape
for (i in 1:length(Ontogeny_list)){ 
  temp_stage <- Ontogeny_list[i]
  temp_grep <- grep(temp_stage, lateral_covariates$Ontogeny)
  temp_lm <- procD.lm(Croc_Lateral_GPA$coords[,,temp_grep] ~ lateral_covariates$simple_ecomorph[temp_grep])
  
  Croc_Lateral_Mean_shape_LMs[[temp_stage]] <- temp_lm
}
#
lapply(Croc_Lateral_Mean_shape_LMs, "summary") #Summary of differences in mean shape among ecomorphs at each individual stage class
p.adjust(c(0.108,0.002,0.018,0.003,0.001,0.001)) #Correction for multiple comparisons
###

### Test for differences in ontogenetic trajectories of ecomorphs ###
# Whole ontogeny, across ecomorphs
Croc_Lateral_Ecomorph_Allometry_LM <- procD.lm(Croc_Lateral_GPA$coords ~ log(Croc_Lateral_GPA$Csize) * lateral_covariates$simple_ecomorph)
summary(Croc_Lateral_Ecomorph_Allometry_LM)
summary(pairwise(Croc_Lateral_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph, covariate = log(Croc_Lateral_GPA$Csize)), test.type = "VC")
Croc_Lateral_CACs <- plotAllometry(Croc_Lateral_Ecomorph_Allometry_LM, log(Croc_Lateral_GPA$Csize), method="CAC")$CAC #calculation of CAC
#
# Whole ontogeny, across species
Croc_Lateral_Species_Allometry_LM <- procD.lm(Croc_Lateral_GPA$coords ~ log(Croc_Lateral_GPA$Csize) * lateral_covariates$species)
summary(Croc_Lateral_Species_Allometry_LM)
summary(pairwise(Croc_Lateral_Species_Allometry_LM, groups = lateral_covariates$species, covariate = log(Croc_Lateral_GPA$Csize)), test.type = "VC")
#
# Embryonic development only, across ecomorphs
Embryos_Croc_Lateral_Ecomorph_Allometry_LM <- procD.lm(Croc_Lateral_GPA$coords[,,Embryos_list] ~ log(Croc_Lateral_GPA$Csize[Embryos_list]) * lateral_covariates$simple_ecomorph[Embryos_list])
summary(Embryos_Croc_Lateral_Ecomorph_Allometry_LM) #P-ANOVA table
summary(pairwise(Embryos_Croc_Lateral_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[Embryos_list])) #Pairwise differences in mean shape
summary(pairwise(Embryos_Croc_Lateral_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[Embryos_list], covariate = log(Croc_Lateral_GPA$Csize[Embryos_list])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#
# Embryonic development only, across species
Embryos_Croc_Lateral_Species_Allometry_LM <- procD.lm(Croc_Lateral_GPA$coords[,,Embryos_list[-29]] ~ log(Croc_Lateral_GPA$Csize[Embryos_list[-29]]) * lateral_covariates$species[Embryos_list[-29]])
summary(Embryos_Croc_Lateral_Species_Allometry_LM) #P-ANOVA table
# must exclude the Gavilis embryos from pairwise comparisons, as with only a single individual the comparisons are impossible otherwise
summary(pairwise(Embryos_Croc_Lateral_Species_Allometry_LM, groups = lateral_covariates$species[Embryos_list[-29]])) #Pairwise differences in mean shape
summary(pairwise(Embryos_Croc_Lateral_Species_Allometry_LM, groups = lateral_covariates$species[Embryos_list[-29]], covariate = log(Croc_Lateral_GPA$Csize[Embryos_list[-29]])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#
# Post-hatching growth only, across ecomorphs
Posthatching_Croc_Lateral_Ecomorph_Allometry_LM <- procD.lm(Croc_Lateral_GPA$coords[,,-Embryos_list] ~ log(Croc_Lateral_GPA$Csize[-Embryos_list]) * lateral_covariates$simple_ecomorph[-Embryos_list])
summary(Posthatching_Croc_Lateral_Ecomorph_Allometry_LM) #P-ANOVA table
summary(pairwise(Posthatching_Croc_Lateral_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[-Embryos_list])) #Pairwise differences in mean shape
summary(pairwise(Posthatching_Croc_Lateral_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[-Embryos_list], covariate = log(Croc_Lateral_GPA$Csize[-Embryos_list])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#
# Post-hatching growth only, across species
Posthatching_Croc_Lateral_Species_Allometry_LM <- procD.lm(Croc_Lateral_GPA$coords[,,-Embryos_list] ~ log(Croc_Lateral_GPA$Csize[-Embryos_list]) * lateral_covariates$species[-Embryos_list])
summary(Posthatching_Croc_Lateral_Species_Allometry_LM) #P-ANOVA table
summary(pairwise(Posthatching_Croc_Lateral_Species_Allometry_LM, groups = lateral_covariates$species[-Embryos_list])) #Pairwise differences in mean shape
summary(pairwise(Posthatching_Croc_Lateral_Species_Allometry_LM, groups = lateral_covariates$species[-Embryos_list], covariate = log(Croc_Lateral_GPA$Csize[-Embryos_list])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#

#correction of p values for pairwise comparisons
p.adjust(c(0.001,0.001,0.002,0.001,0.001,0.021,0.001,0.001,0.014))
p.adjust(c(0.035,0.016,0.001,0.637,0.203,0.003,0.171,0.060,0.024))
p.adjust(c(0.154,0.056,0.001,0.826,0.366,0.003,0.09,0.02,0.02)) #Correction for multiple comparisons among ecomorphs
###


###Calculate slopes for each trajectory#
Croc_Lateral_Slopes <- list()
# for loop to iterate through each ontogenetic class and perform P-ANOVA for ecomorph differences in mean shape
for(j in 1:3){
  temp_group <- c("Complete", "Embryonic", "Post-hatching")[j]
  temp_CAC <- list(Croc_Lateral_CACs, Croc_Lateral_CACs[Embryos_list], Croc_Lateral_CACs[-Embryos_list])[[j]]
  temp_CS <- list(Croc_Lateral_GPA$Csize, Croc_Lateral_GPA$Csize[Embryos_list], Croc_Lateral_GPA$Csize[-Embryos_list])[[j]]
  temp_cov <- list(lateral_covariates, lateral_covariates[Embryos_list,], lateral_covariates[-Embryos_list,])[[j]]
  slope_matrix <- matrix(data=NA,nrow=3,ncol=2,dimnames = list(c(unique(lateral_covariates$simple_ecomorph)),c("Intercept","Slope")))
  for (i in 1:length(unique(lateral_covariates$simple_ecomorph))){ 
    temp_eco <- unique(temp_cov$simple_ecomorph)[i]
    temp_grep <- grep(temp_eco, temp_cov$simple_ecomorph)
    temp_lm <- lm(temp_CAC[temp_grep] ~ log(temp_CS[temp_grep]))
    slope_matrix[temp_eco,] <- temp_lm$coefficients
  }
  Croc_Lateral_Slopes[[temp_group]] <- slope_matrix
}
Croc_Lateral_Slopes
###