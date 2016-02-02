############################################################################
#LEMUR BODY MASS VERSUS MAX INGESTED FOOD SIZE
############################################################################
#load libraries
library(xlsx);library(ape);library(phytools);library(coda);library(diversitree);library(geiger);library(nlme);library(phylolm)
#######################################################################
#set working directory
getwd()
#set your working directory. yours will be different from mine.
setwd("~/Dropbox/Canareae_project/functional morphology/analyses for supplement/LEMURS/REGRESSION")
#####################################################################
#read in character data for those lemurs that are frugivorous

topredict<-read.csv("Body_mass_for_Vb_estimates.csv") # these are the lemur species for which we need to predict Vb
vb<-read.csv("empirical_Vb_measurements.csv") #these are the lemur species that have empirical Vb measurements. We will use these data to predict the Vb of the other species
#Log transform data
vb$mass<-log(vb[,2]) #this log transforms body mass
vb$vb<-log(vb[,3]) #this log transforms the empirical Vb measurements
topredict$mass<-log(topredict[,2]) #this log transforms body mass of species for which we're predicting Vb
######################################################################
#run a linear regression of lg Vb ~ log body mass to establish the predictive value
trait<-lm(vb~mass, data=vb)
summary(trait)
coefficients(trait)
plot(vb~mass, data=vb, pch=16)
abline(trait, col="black")
#########################################################################
#OK, so we've established a strong linear relationship between body mass and Vb. Now, let's predict the Vb of the remaining frugivorous lemurs in our phylogeny with their body mass.
plot(vb~mass, data=vb, ylim=c(-0.6,2.5), xlim=c(-3,4), col="black", pch=16)
abline(trait, col="black")
###make new dataframe 
output=matrix(nrow=length(topredict$mass), ncol=3)
for (i in 1:length(topredict$mass))
{
	newdata=data.frame(mass=topredict$mass[i])
predict(trait, newdata,interval="predict")->output[i,]
	}
	
output #these are your predicted vb values and their standard deviations
points(topredict$mass,output[,1], col= "red", pch=16) #this maps your predicted Vb values onto your regression line

#add lines to your regression showing the log average seed diameter of putatively orphaned lineages, and the 95% confidence interval surrounding Canarium fruit diameter (black lines)

abline(h=0.929, col="black", lty=4 ) 
abline(h=(0.783), col="black" )
abline(h=(1.078), col="black")
abline(h=(0.811), col="grey" )
abline(h=(1.31), col="grey" )
abline(h=(0.99), col="grey" )
abline(h=(1.16), col="grey" )
abline(h=(1.93), col="grey" )
abline(h=(0.405), col="grey" )
