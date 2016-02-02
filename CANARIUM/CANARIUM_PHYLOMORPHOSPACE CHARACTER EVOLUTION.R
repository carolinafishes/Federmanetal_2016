###################################################
#PHYLOMORPHOSPACE & CHARACTER EVOLUTION, CANARIEAE FRUIT
###################################################
#load libraries
library(ape);library(phytools);library(coda);library(diversitree);library(geiger)
###################################################
#set working directory, this will be different for each user.
getwd()
setwd("~/Dropbox/Canareae_project/functional morphology/analyses for supplement/CANARIUM")
###################################################
#MALAGASY CANARIUM PHYLOMORPHOSPACE
###################################################
#read in character data
malagasy<-read.csv("malagasy_canarium_fruit.csv")
malagasy.data<-data.frame(malagasy[,2:3])
malagasy.data$length<-log(malagasy.data[,1]/10)
malagasy.data$width<-log(malagasy.data[,2]/10)
malagasy.data<-malagasy.data[,3:4]
rownames(malagasy.data)<-malagasy[,1]
#read in tree
canarieae<-read.tree("Canarieae_ingroup.phy")
##################################################
#check overlap of phylogenetic and trait datasets
name.check(canarieae, malagasy.data)->phylOverlap
phylOverlap$tree_not_data
malagasy_canarium<-drop.tip(canarieae, phylOverlap$tree_not_data)
plot(malagasy_canarium)
name.check(malagasy_canarium, malagasy.data)
######################################################
#######################################################################
#time color coded morphospace
x<-as.matrix(malagasy.data)[,1]
AA<-contMap(malagasy_canarium, x)
 H<-nodeHeights(malagasy_canarium)
 h<-max(H)
 
 
 
  AA$cols[]<-rainbow(1001,start=0.7,end=4/6)
   for(i in 1:nrow(H)) names(AA$tree$maps[[i]])<- round((H[i,1]+cumsum(AA$tree$maps[[i]]))/h*1000)

 plot(AA,legend=FALSE)
 # check to verify that temporal information is correct
 #plot(AA,legend=FALSE)
 
 mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 


 phylomorphospace(AA$tree,malagasy.data,colors=AA$cols,lwd=3, node.by.map=TRUE,ylim=c(0.74,1.38), xlim=c(1.0,1.65))
 
#############################################################################################
#ANCESTRAL STATE RECONSTRUCTION OF CANARIEAE FRUIT WIDTH
#############################################################################################
# import and log-transformed fruit width for all of the Canarieae
########
fruit<-read.csv("Canarieae_fruit.csv")
names(fruit)
fruit.data<-data.frame(fruit[,3])
fruit.data$width<-log(fruit.data[,1]/10)
rownames(fruit.data)<-fruit[,1]

#########################

#now, re-prune the Canarieae tree to include all of the Canarieae from the dataset
name.check(canarieae, fruit.data)->phylOverlap phylOverlap$tree_not_data
canarieae<-drop.tip(phyl, phylOverlap$tree_not_data)
plot(canTree)
name.check(canTree, fruit.data)



############
#pull out the character of interest

width <-data.frame(fruit.data[,2])
rownames(width)<-fruit[,1]
X<-as.matrix(width, row.names=1)[,1]
head(X)

############
#do the analysis
obj<-contMap(canTree, X, res=100, plot=T, method="anc.ML", lwd=0.8, fsize=0.55, outline=F)
