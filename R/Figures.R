#### Figures ####

### Morphospace Figures ###
PCA <- Croc_Lateral_PCA #change to produce figures for other datasets
PCA$importance <- PCA$d/sum(PCA$d)
PV <- lateral_covariates
PC_1 <- 1
PC_2 <- 2

Xlim<-c(floor(min(PCA$x[,PC_1])*10)/10,ceiling(max(PCA$x[,PC_1])*10)/10)
Ylim<-c(floor(min(PCA$x[,PC_2])*10)/10,ceiling(max(PCA$x[,PC_2])*10)/10)

plot(0, 0, type = "n",
     xlim = Xlim,
     ylim = Ylim,
     xlab = paste("Principal Component ", PC_1, " (", round(100*PCA$importance[PC_1], digits = 1), "%)", sep = ""),
     ylab = paste("Principal Component ", PC_2, " (", round(100*PCA$importance[PC_2], digits = 1), "%)", sep = ""),
     axes = FALSE,
     frame.plot = FALSE,
     asp=F)

axis(1, round(seq(Xlim[1],Xlim[2],by=0.1),1), pos=Ylim[1])
axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1], las = 1)
clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
abline(h=0, lty=3)
abline(v=0, lty=3)

clip(-0.4,0.4,-0.4,0.4)

points(PCA$x[,PC_1],PCA$x[,PC_2],
       cex = PV$Size,
       pch = PV$Shape,
       bg = alpha(PV$Color, 0.75), asp=F)

###

### Allometric Figures ###
GPA <- Croc_Lateral_GPA
CACs <- Croc_Lateral_CACs #change to produce figures for other datasets
PV <- lateral_covariates

Xlim<-c(0.0,ceiling(max(log(GPA$Csize))))
Ylim<-c(floor(min(CACs)*10)/10,ceiling(max(CACs)*10)/10)

plot(0, 0, type = "n",
     xlim = Xlim,
     ylim = Ylim,
     xlab = paste("Centroid Size (log scale)"),
     ylab = paste("Common Allometric Component"),
     axes = FALSE,
     frame.plot = FALSE,
     asp=F)

axis(1, round(seq(Xlim[1],Xlim[2],by=0.5),1), pos=Ylim[1])
axis(2, round(seq(Ylim[1],Ylim[2],by=0.1),1), pos=Xlim[1], las = 1)
clip(Xlim[1],Xlim[2],Ylim[1],Ylim[2])
abline(h=0, lty=3)

clip(0,5,-0.4,0.4)

points(log(GPA$Csize), CACs,
       cex = PV$Size,
       pch = PV$Shape,
       bg = alpha(PV$Color, 0.75), asp=F)
###