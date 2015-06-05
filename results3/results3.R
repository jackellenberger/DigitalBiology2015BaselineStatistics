r55_Dehydrons<-read.csv(file="radius55/perDehydron_r55.csv",header=TRUE)
r55_HBonds<-read.csv(file="radius55/perHBond_r55.csv",header=TRUE)
r55_PerFile<-read.csv(file="radius55/perFile_r55.csv",header=TRUE)
r55_Model<-read.csv(file="radius55/ModelDehydrons_r55.csv",header=TRUE)

r60_Dehydrons<-read.csv(file="radius60/perDehydron_r60.csv",header=TRUE)
r60_HBonds<-read.csv(file="radius60/u_perHBond_r60.csv",header=TRUE)
r60_PerFile<-read.csv(file="radius60/u_perFile_r60.csv",header=TRUE)
r60_Model<-read.csv(file="radius60/ModelDehydrons_r60.csv",header=TRUE)

r65_Dehydrons<-read.csv(file="radius65/perDehydron_r65.csv",header=TRUE)
r65_HBonds<-read.csv(file="radius65/perHBond_r65.csv",header=TRUE)
r65_PerFile<-read.csv(file="radius65/perFile_r65.csv",header=TRUE)
r65_Model<-read.csv(file="radius65/ModelDehydrons_r65.csv",header=TRUE)

r70_Dehydrons<-read.csv(file="radius70/perDehydron_r70.csv",header=TRUE)
r70_HBonds<-read.csv(file="radius70/u_perHBond_r70.csv",header=TRUE)
r70_PerFile<-read.csv(file="radius70/u_perFile_r70.csv",header=TRUE)
r70_Model<-read.csv(file="radius70/ModelDehydrons_r70.csv",header=TRUE)

r75_Dehydrons<-read.csv(file="radius75/perDehydron_r75.csv",header=TRUE)
r75_HBonds<-read.csv(file="radius75/perHBond_r75.csv",header=TRUE)
r75_PerFile<-read.csv(file="radius75/perFile_r75.csv",header=TRUE)
r75_Model<-read.csv(file="radius75/ModelDehydrons_r75.csv",header=TRUE)

ddatoms<-read.csv("../results4/r_nrpdb_1_modeledAtomsInDD_output.csv",header=TRUE)
pfile<-read.csv("../results4/r_nrpdb_1_perFile_output.csv",header=TRUE)

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


##### figure 5: expected vs actual dehydron count #####
## 6.5A ##########
hb1_65<-remove_outliers(r65_PerFile$HBless1SD)
dc_65<-remove_outliers(r65_PerFile$dehydronCount)
plot(hb1_65,dc_65,
     na.rm=TRUE,
     main="Count of Actual vs Predicted Dehydrons",
     pch=20, cex=.4,
     xlab="# Predicted Dehydrons (1 S.D. Below Mean # Hydrogen Bonds)",
     ylab="# Actual Dehydrons"
     )
mtext("Desolvation Domain Radius = 6.5A")

## 6.5A ##########
hb1_65<-remove_outliers(r65_PerFile$HBless1SD)
dc_65<-remove_outliers(r65_PerFile$dehydronCount)
plot(hb1_65,dc_65,
     na.rm=TRUE,
     main="Count of Actual vs Predicted Dehydrons",
     pch=20, cex=.4,
     xlab="# Predicted Dehydrons (1 S.D. Below Mean # Hydrogen Bonds)",
     ylab="# Actual Dehydrons"
     )
mtext("Desolvation Domain Radius = 6.5A, y=x")
abline(0,1)

## 6.0A ##########
hb1_60<-remove_outliers(r60_PerFile$HBless1SD)
dc_60<-remove_outliers(r60_PerFile$dehydronCount)
plot(hb1_60,dc_60,
     na.rm=TRUE,
     main="Count of Actual vs Predicted Dehydrons",
     pch=20, cex=.4,
     xlab="# Predicted Dehydrons (1 S.D. Below Mean # Hydrogen Bonds)",
     ylab="# Actual Dehydrons"
     )
mtext("Desolvation Domain Radius = 6.0A, y=x")
abline(0,1)

## 7.0A ##########
hb1_70<-remove_outliers(r70_PerFile$HBless1SD)
dc_70<-remove_outliers(r70_PerFile$dehydronCount)
plot(hb1_70,dc_70,
     na.rm=TRUE,
     main="Count of Actual vs Predicted Dehydrons",
     pch=20, cex=.4,
     xlab="# Predicted Dehydrons (1 S.D. Below Mean # Hydrogen Bonds)",
     ylab="# Actual Dehydrons"
     )
mtext("Desolvation Domain Radius = 7.0A, y=x")
abline(0,1)

##### figure 6 Left #####
## 5.5A ##########
f_Addd_55<-r55_Dehydrons$uniqueAtomsInDDD[r55_Dehydrons$uniqueAtomsInDDD<300]
f_Mddd_55<-r55_Model$uniqueAtomsInDDD[r55_Model$uniqueAtomsInDDD<300]
Addd_55<-head(f_Addd_55,5000)
Mddd_55<-head(f_Mddd_55,5000)
plot(Addd_55,Mddd_55,
     xlim=c(0,300),ylim=c(0,350),
     pch=20,cex=.1,
     main="Number of Atoms in Desolvation Domain, Observed vs modeled",
     xlab="Number of Atoms in Observed Dehydron DD",
     ylab="Number of Atoms in modeled Dehydron DD"
)
mtext("Desolvation Domain Radius = 5.5A, with y=x")
abline(0,1)
points(mean(Addd_55),mean(Mddd_55),pch=20,col="red")


## 6.0A ##########
Addd_60<-head(r60_Dehydrons$uniqueAtomsInDDD,5000)
Mddd_60<-head(r60_Model$uniqueAtomsInDDD,5000)
plot(Addd_60,Mddd_60,
     xlim=c(0,300),ylim=c(0,350),
     pch=20,cex=.1,
     main="Number of Atoms in Desolvation Domain, Observed vs modeled",
     xlab="Number of Atoms in Observed Dehydron DD",
     ylab="Number of Atoms in modeled Dehydron DD"
)
mtext("Desolvation Domain Radius = 6.0A, with y=x")
abline(0,1)
points(mean(Addd_60),mean(Mddd_60),pch=20,col="red")


## 6.5A ##########
Addd_65<-head(r65_Dehydrons$uniqueAtomsInDDD,5000)
Mddd_65<-head(r65_Model$uniqueAtomsInDDD,5000)
plot(Addd_65,Mddd_65,
     xlim=c(0,300),ylim=c(0,350),
     pch=20,cex=.1,
     main="Number of Atoms in Desolvation Domain, Observed vs modeled",
     xlab="Number of Atoms in Observed Dehydron DD",
     ylab="Number of Atoms in modeled Dehydron DD"
)
mtext("Desolvation Domain Radius = 6.5A, with y=x")
abline(0,1)
points(mean(Addd_65),mean(Mddd_65),pch=20,col="red")


## 7.0A ##########
Addd_70<-head(r70_Dehydrons$uniqueAtomsInDDD,5000)
Mddd_70<-head(r70_Model$uniqueAtomsInDDD,5000)
plot(Addd_70,Mddd_70,
     xlim=c(0,300),ylim=c(0,350),
     pch=20,cex=.1,
     main="Number of Atoms in Desolvation Domain, Observed vs modeled",
     xlab="Number of Atoms in Observed Dehydron DD",
     ylab="Number of Atoms in modeled Dehydron DD"
)
mtext("Desolvation Domain Radius = 7.0A, with y=x")
abline(0,1)
points(mean(Addd_70),mean(Mddd_70),pch=20,col="red")


## 7.5A ##########
Addd_75<-head(r75_Dehydrons$uniqueAtomsInDDD,2800)
Mddd_75<-head(r75_Model$uniqueAtomsInDDD,2800)
plot(Addd_75,Mddd_75,
     xlim=c(0,300),ylim=c(0,350),
     pch=20,cex=.1,
     main="Number of Atoms in Desolvation Domain, Observed vs modeled",
     xlab="Number of Atoms in Observed Dehydron DD",
     ylab="Number of Atoms in modeled Dehydron DD"
)
mtext("Desolvation Domain Radius = 7.5A, with y=x")
abline(0,1)
points(mean(Addd_75),mean(Mddd_75),pch=20,col="red")

plot(Addd_55,Mddd_55,pch=20,cex=.06,col="violet",
     xlab="Atoms in Observed Dehydron Desolvation Domain",
     ylab="Atoms in Modeled Dehydron Desolvation Domain",
     xlim=c(50,250))
points(Addd_60,Mddd_60,pch=20,cex=.06,col="blue")
points(Addd_65,Mddd_65,pch=20,cex=.06,col="green")
points(Addd_70,Mddd_70,pch=20,cex=.06,col="orange")
points(Addd_75,Mddd_75,pch=20,cex=.06,col="red")
points(mean(Addd_55),mean(Mddd_55),pch=20,cex=2,col="darkviolet")
points(mean(Addd_60),mean(Mddd_60),pch=20,cex=2,col="darkblue")
points(mean(Addd_65),mean(Mddd_65),pch=20,cex=2,col="darkgreen")
points(mean(Addd_70),mean(Mddd_70),pch=20,cex=2,col="darkorange")
points(mean(Addd_75),mean(Mddd_75),pch=20,cex=2,col="darkred")
legend("bottomright", c("5.5A", "6.0A","6.5A","7.0A","7.5A"), col=c("darkviolet","darkblue","darkgreen","darkorange","darkred"), lwd=4)

plot.new()
##### Figure 6 Right#####
## 6.5A #####
#f_Mddd_65 <- Mddd_65[Mddd_65<250]
breaks=c(0,25,50,75,100,125,150,175,200,225,250,275,300,325)
h_Addd_65<-r65_Dehydrons$uniqueAtomsInDDD
h_Mddd_65<-r65_Model$uniqueAtomsInDDD
plot_Addd_65<-hist(h_Addd_65,breaks=1000)
plot_Mddd_65<-hist(h_Mddd_65[h_Mddd_65<275],breaks=1000)
plot(plot_Mddd_65, border=rgb(1,0,0,1/4), col=rgb(1,0,0,1/4), xlim=c(50,250),ylim=c(0,830),
	main="Atoms in Observed (blue) and Modeled (red) Desolvation Domains",
	xlab="Number of Atoms in Desolvation Domain")
plot(plot_Addd_65, border=rgb(0,0,1,1/4), col=rgb(0,0,1,1/4),xlim=c(50,250),ylim=c(0,830), add=T)
legend("topright", c("Observed", "Modeled"), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), lwd=4)
mtext("Desolvation Domain Radius = 6.5A, per dehydron")


## 6.0A #####
f_Mddd_60 <- Mddd_60[Mddd_60<250]
breaks=c(0,25,50,75,100,125,150,175,200,225,250,275,300,325)
plot_Addd_60<-hist(Addd_60,breaks=1000)
plot_Mddd_60<-hist(f_Mddd_60,breaks=1000)
plot(plot_Mddd_60, border=rgb(1,0,0,1/4), col=rgb(1,0,0,1/4), xlim=c(50,250),ylim=c(0,150),
	main="Atoms in Observed (blue) and Modeled (red) Desolvation Domains",
	xlab="Number of Atoms in Desolvation Domain")
plot(plot_Addd_60, border=rgb(0,0,1,1/4), col=rgb(0,0,1,1/4),xlim=c(50,250),ylim=c(0,150), add=T)
legend("topright", c("Observed", "Modeled"), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), lwd=4)
mtext("Desolvation Domain Radius = 6.0A")


## 7.0A #####
f_Mddd_70 <- Mddd_70[Mddd_70<250]
breaks=c(0,25,50,75,100,125,150,175,200,225,250,275,300,325)
plot_Addd_70<-hist(Addd_70,breaks=1000)
plot_Mddd_70<-hist(f_Mddd_70,breaks=1000)
plot(plot_Mddd_70, border=rgb(1,0,0,1/4), col=rgb(1,0,0,1/4), xlim=c(50,250),ylim=c(0,120),
	main="Atoms in Observed (blue) and Modeled (red) Desolvation Domains",
	xlab="Number of Atoms in Desolvation Domain")
plot(plot_Addd_70, border=rgb(0,0,1,1/4), col=rgb(0,0,1,1/4),xlim=c(50,250),ylim=c(0,120), add=T)
legend("topright", c("Observed", "Modeled"), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), lwd=4)
mtext("Desolvation Domain Radius = 7.0A")

########## VISHAL's PLOTS #####

breaks2=c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000)
modeledAtomHist<-hist(ddatoms$modeledUniqueAtomsInDDD[ddatoms$modeledUniqueAtomsInDDD<6000],breaks=breaks2)
actualAtomHist<-hist(ddatoms$actualUniqueAtomsInDDD,breaks=breaks2)
plot(modeledAtomHist,border=rgb(1,0,1,1/4),col=rgb(1,0,1,1/4),xlim=c(0,6000),ylim=c(0,300),main="Actual (blue) vs Modeled (red) # Unique Atoms in DDD",xlab="# Anique atoms in DDD")
plot(actualAtomHist,border=rgb(0,0,1,1/4),col=rgb(0,0,1,1/4),xlim=c(0,6000),add=T)
mtext("Desolvation Domain Radius = 6.5A, Per Protein")
legend("topright", c("Observed", "Modeled"), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), lwd=4)

names(ddatoms)
diffsM<-remove_outliers(abs(ddatoms$numNonUniqueAtomsInDDD - ddatoms$modeledUniqueAtomsInDDD))
mean(diffsM,na.rm = TRUE)
median(diffsM,na.rm=TRUE)

names(pfile)
diffsA<-remove_outliers(abs(pfile$totalAtomsInDDD - pfile$totalUniqueAtomsInDDD))
mean(diffsA,na.rm = TRUE)
median(diffsA,na.rm = TRUE)


wrapPercentM <- ddatoms$modeledUniqueAtomsInDDD / pfile$totalAtomsInPDB
mean(wrapPercentM)
wrapPercentA <- pfile$totalUniqueAtomsInDDD / pfile$totalAtomsInPDB
mean(wrapPercent)

