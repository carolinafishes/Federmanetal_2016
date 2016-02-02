#######################################################################
#LEMUR MORPHOSPACE
######################################################################
#load libraries
library(ape);library(phytools);library(coda);library(diversitree);library(geiger);library(nlme);library(phylolm);library(PBSmapping);library(siar);library(spatstat)
#######################################################################
#set working directory
getwd()
#set your working directory. This will be different for each user.
setwd("~/Dropbox/Canareae_project/functional morphology/analyses for supplement/LEMURS/PHYLOMORPHOSPACE")
#####################################################################
#read in character data and log-transform
traits<-read.csv("phylomorphospace_data.csv")
trait.data<-data.frame(traits) 
rownames(trait.data)<-traits[,1] #make species rownames
names(trait.data)
trait.data<-trait.data[,5:14] #define the trait data
row.names(trait.data)
size<-c(traits[,4]) #pull out body mass
log(size)->lsize
g<-unique(traits[,1])
names(lsize)<-g
log(trait.data)->ltrait.data
ltrait<-as.matrix(ltrait.data)
##################################################################
#read in tree and prune out species that aren't in the trait dataset (i.e. the outgroups)
lemur_tree<-read.nexus ("Final_lemur")
name.check(lemur_tree, lsize)->lemurOverlap
lemurOverlap$tree_not_data
lemurTree<-drop.tip(lemur_tree, lemurOverlap$tree_not_data)
name.check(lemurTree, lsize)
##################################################################
##phylogenetic regressions to get the residuals
phyl.resid(lemurTree,lsize,ltrait,method="BM")->resids
as.matrix(resids$resid)->X
head(X)
result<-phyl.pca(lemurTree, X,method="lambda",mode="corr")
result$S->fruitomorphospaceAll
summary(result)
result$S
#######################################################################
#time color colded phylomorphospace

AA<-contMap(lemurTree, fruitomorphospaceAll[,1])
 H<-nodeHeights(lemurTree)
 h<-max(H)
 
 
 
  AA$cols[]<-rainbow(1001,start=0.7,end=4/6)
   for(i in 1:nrow(H)) names(AA$tree$maps[[i]])<- round((H[i,1]+cumsum(AA$tree$maps[[i]]))/h*1000)

 plot(AA,legend=FALSE)
 # check to verify that temporal information is correct
 #plot(AA,legend=FALSE)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 

phylomorphospace(AA$tree,fruitomorphospaceAll[,c(1,2)],colors=AA$cols,lwd=3, node.by.map=TRUE,xlab="trait 1",ylab="trait 2 ")
 



######################################################################
#check overlap of fruitogenetic and trait datasets
name.check(fruit, trait.data)->fruitOverlap
fruitOverlap$tree_not_data
lemurTree<-drop.tip(fruit, fruitOverlap$tree_not_data)
plot(lemurTree)
write.tree(lemurTree, "lemurTree")
##########################################################################################
#STOCHASTIC CHARACTER MAPPING
#########################################################################################
names(traits)
traits[,2]
feeding <-data.frame(traits[,2])
rownames(feeding)<-traits[,1]
#
name.check(lemur_tree, feeding)->fruitOverlap
fruitOverlap$tree_not_data
lemurTree<-drop.tip(lemur_tree, lemurOverlap$tree_not_data)
#make sure that tip names and trait names are in the same order
data <- feeding[match(lemurTree$tip.label,rownames(feeding)),]
X<-as.matrix(feeding, row.names=1)[,1]
mtreez<-make.simmap(lemurTree, X, nsim=1000)
colors<-setNames(c("blue", "coral", "green", "purple", "gold1"), sort(unique(getStates(mtreez[1], "tips"))))
plotSimmap(mtreez[[1]], colors=colors, pts=F, ftype="off")
plot(lemurTree, mar=c(5.1,0.2,0.2,0.2))
obj<-describe.simmap(mtreez)
nodelabels(pie=obj$ace, piecol=colors, cex=0.6)
##########################################################################################



