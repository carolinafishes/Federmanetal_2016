#####################################################################################################
#CALCULATE THE AREA OF PHYLOMORPHOSPACE LOST WITH THE EXTINCTION OF LARGE-BODIED DISPERSERS
######################################################################################################
#set up your workspace, close all currently open windows, set the working directory and load libraries 
graphics.off()
getwd()
setwd("~/Dropbox/Canareae_project/functional morphology/analyses for supplement/LEMURS/PHYLOMORPHOSPACE")
library(xlsx);library(ape);library(phytools);library(coda);library(diversitree);library(geiger);library(nlme);library(MCMCglmm);library(phylolm);library(PBSmapping);library(siar);library(spatstat);library("sp");library("rgdal"); library(rgeos)
######################################################################################################
# read your phylogenetic PCA of the residuals for probable dispersers (here, only PC axes 1 and 2)
# the column names have to be exactly, "group", "x", "y"(here, grouping is by extinct: group 2, and extant: group 1)

fruit<-read.csv("frugivores_pca.csv")
ext<-read.csv("extant.csv")

# split the extinct and extant taxa data based on group
spx <- split(fruit$x,fruit$Group)
spy <- split(fruit$y,fruit$Group)

#now draw convex hulls around your data, subsetting it by extinct, extant, and everything together
extant<-convexhull(spx[[1]],spy[[1]])
extant <- c(extant, extant[1])
all<-convexhull(fruit$x,fruit$y)
all <- c(all, all[1])
extinct<-convexhull(spx[[2]], spy[[2]] )
extinct <- c(extinct, extinct[1])
#now plot them so you can see

plot(all$xcoords, all$ycoords, col="blue")
lines(all$xcoords, all$ycoords, col="blue")
points(extinct$xcoords, extinct$ycoords, col="gold1")
lines(extinct$xcoords, extinct$ycoords, col="gold1")
points(extant$xcoords, extant$ycoords, col="green")
lines(extant$xcoords, extant$ycoords, col="green")

#subtract the area of extant taxa from the area of all taxa to get the area lost
lostArea<-all$TA-extant$TA # lost area is 248.9742
lostArea

#now divide the lostArea by the total area and * 100 to get the percentage of morphospace lost

pct<-(lostArea/all$TA)*100
pct



