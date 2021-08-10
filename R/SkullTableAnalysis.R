### Skull Table Analysis ###


### Perform GPA alignment for all three datasets ###
SkullTable_20_GPA <- gpagen(Croc_SkullTable_SubCurve_20,
                            curves=Curves_20pt_sliders,
                            ProcD=TRUE)
###

### Perform PCA for all three datasets ###
SkullTable_20_PCA <- gm.prcomp(SkullTable_20_GPA$coords) #skull table only dataset PCA
###

### Test correclation between PCs and centroid size ###
cor(SkullTable_20_PCA$x,log(Croc_Lateral_GPA$Csize))

p.cor <- NULL
for (i in 1:ncol(SkullTable_20_PCA$x)){
  p.cor <- rbind(p.cor,cor.test(log(Croc_Lateral_GPA$Csize),SkullTable_20_PCA$x[,i])$p.value)
}
###

### Test for mean shape differences among ecomorphs at each ontogenetic stage ###
SkullTable_20_Mean_shape_LMs <- list()
# for loop to iterate through each ontogenetic class and perform P-ANOVA for ecomorph differences in mean shape
for (i in 1:length(Ontogeny_list)){ 
  temp_stage <- Ontogeny_list[i]
  temp_grep <- grep(temp_stage, lateral_covariates$Ontogeny)
  temp_lm <- procD.lm(SkullTable_20_GPA$coords[,,temp_grep] ~ lateral_covariates$simple_ecomorph[temp_grep])
  
  SkullTable_20_Mean_shape_LMs[[temp_stage]] <- temp_lm
}
#
lapply(SkullTable_20_Mean_shape_LMs, "summary") #Summary of differences in mean shape among ecomorphs at each individual stage class
p.adjust(c(0.002,0.175,0.72,0.102,0.02,0.089)) #Correction for multiple comparisons
###

### Test for differences in ontogenetic trajectories of ecomorphs & species ###
# Whole ontogeny, across ecomorphs
SkullTable_Ecomorph_Allometry_LM <- procD.lm(SkullTable_20_GPA$coords ~ log(Croc_Lateral_GPA$Csize) * lateral_covariates$simple_ecomorph)
summary(SkullTable_Ecomorph_Allometry_LM)
summary(pairwise(SkullTable_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph, covariate = log(Croc_Lateral_GPA$Csize)), test.type = "VC")
SkullTable_20_CACs <- plotAllometry(SkullTable_Ecomorph_Allometry_LM, log(Croc_Lateral_GPA$Csize), method="CAC")$CAC #calculation of CAC
#
# Whole ontogeny, across species
SkullTable_Species_Allometry_LM <- procD.lm(SkullTable_20_GPA$coords ~ log(Croc_Lateral_GPA$Csize) * lateral_covariates$species)
summary(SkullTable_Species_Allometry_LM)
summary(pairwise(SkullTable_Species_Allometry_LM, groups = lateral_covariates$species, covariate = log(Croc_Lateral_GPA$Csize)), test.type = "VC")
#
# Embryonic development only, across ecomorphs
Embryos_SkullTable_Ecomorph_Allometry_LM <- procD.lm(SkullTable_20_GPA$coords[,,Embryos_list] ~ log(Croc_Lateral_GPA$Csize[Embryos_list]) * lateral_covariates$simple_ecomorph[Embryos_list])
summary(Embryos_SkullTable_Ecomorph_Allometry_LM) #P-ANOVA table
summary(pairwise(Embryos_SkullTable_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[Embryos_list])) #Pairwise differences in mean shape
summary(pairwise(Embryos_SkullTable_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[Embryos_list], covariate = log(Croc_Lateral_GPA$Csize[Embryos_list])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#
# Embryonic development only, across species
Embryos_SkullTable_Species_Allometry_LM <- procD.lm(SkullTable_20_GPA$coords[,,Embryos_list[-29]] ~ log(Croc_Lateral_GPA$Csize[Embryos_list[-29]]) * lateral_covariates$species[Embryos_list[-29]])
summary(Embryos_SkullTable_Species_Allometry_LM) #P-ANOVA table
# must exclude the Gavilis embryos from pairwise comparisons, as with only a single individual the comparisons are impossible otherwise
summary(pairwise(Embryos_SkullTable_Species_Allometry_LM, groups = lateral_covariates$species[Embryos_list[-29]])) #Pairwise differences in mean shape
summary(pairwise(Embryos_SkullTable_Species_Allometry_LM, groups = lateral_covariates$species[Embryos_list[-29]], covariate = log(Croc_Lateral_GPA$Csize[Embryos_list[-29]])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#
# Post-hatching growth only, across ecomorphs
Posthatching_SkullTable_Ecomorph_Allometry_LM <- procD.lm(SkullTable_20_GPA$coords[,,-Embryos_list] ~ log(Croc_Lateral_GPA$Csize[-Embryos_list]) * lateral_covariates$simple_ecomorph[-Embryos_list])
summary(Posthatching_SkullTable_Ecomorph_Allometry_LM) #P-ANOVA table
summary(pairwise(Posthatching_SkullTable_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[-Embryos_list])) #Pairwise differences in mean shape
summary(pairwise(Posthatching_SkullTable_Ecomorph_Allometry_LM, groups = lateral_covariates$simple_ecomorph[-Embryos_list], covariate = log(Croc_Lateral_GPA$Csize[-Embryos_list])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#
# Post-hatching growth only, across species
Posthatching_SkullTable_Species_Allometry_LM <- procD.lm(SkullTable_20_GPA$coords[,,-Embryos_list] ~ log(Croc_Lateral_GPA$Csize[-Embryos_list]) * lateral_covariates$species[-Embryos_list])
summary(Posthatching_SkullTable_Species_Allometry_LM) #P-ANOVA table
summary(pairwise(Posthatching_SkullTable_Species_Allometry_LM, groups = lateral_covariates$species[-Embryos_list])) #Pairwise differences in mean shape
summary(pairwise(Posthatching_SkullTable_Species_Allometry_LM, groups = lateral_covariates$species[-Embryos_list], covariate = log(Croc_Lateral_GPA$Csize[-Embryos_list])), test.type = "VC") #Pairwise differences in slope of ontogenetic trajectories
#

#correction of p values for pairwise comparisons
p.adjust(c(0.092,0.174,0.007,0.249,0.251,0.029,0.150,0.725,0.327))
p.adjust(c(0.365,0.437,0.142,0.302,0.940,0.007,0.466,0.401,0.094)) #Correction for multiple comparisons among ecomorphs
###


###Calculate slopes for each trajectory#
SkullTable_Slopes <- list()
# for loop to iterate through each ontogenetic class and perform P-ANOVA for ecomorph differences in mean shape
for(j in 1:3){
  temp_group <- c("Complete", "Embryonic", "Post-hatching")[j]
  temp_CAC <- list(SkullTable_20_CACs, SkullTable_20_CACs[Embryos_list], SkullTable_20_CACs[-Embryos_list])[[j]]
  temp_CS <- list(Croc_Lateral_GPA$Csize, Croc_Lateral_GPA$Csize[Embryos_list], Croc_Lateral_GPA$Csize[-Embryos_list])[[j]]
  temp_cov <- list(lateral_covariates, lateral_covariates[Embryos_list,], lateral_covariates[-Embryos_list,])[[j]]
  slope_matrix <- matrix(data=NA,nrow=3,ncol=2,dimnames = list(c(unique(lateral_covariates$simple_ecomorph)),c("Intercept","Slope")))
for (i in 1:length(unique(lateral_covariates$simple_ecomorph))){ 
  temp_eco <- unique(temp_cov$simple_ecomorph)[i]
  temp_grep <- grep(temp_eco, temp_cov$simple_ecomorph)
  temp_lm <- lm(temp_CAC[temp_grep] ~ log(temp_CS[temp_grep]))
  slope_matrix[temp_eco,] <- temp_lm$coefficients
}
SkullTable_Slopes[[temp_group]] <- slope_matrix
}
SkullTable_Slopes
###
