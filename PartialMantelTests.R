#Clear the working environment
rm(list=ls())

#I am loading in the packages that are needed for Partial Mantel Tests
library(vegan)#This package contains the vegdist function which is used to make distance matrices from our diet and mycobiome tables
library(ecodist)#This package contains the mantel function which is used to run the partial Mantel tests

#Setting the working directory, this will vary by machine
setwd("D:/01 Documents/School/ND/Macaque Diet/Publications/Diet Mycobiome/Data")

#Loading in the various data files
#This is a file containing sample level descriptions of the fungal communities at a genus resolution for samples from Singapore
SingaporeFungusTable<-read.table("SingaporeFungusGenusCutOffFiltered.txt", quote="",header=T,sep="\t",row.names=1)

#This is a file containing a sample level description of the diet at a genus resolution for samples from Singapore
SingaporeDietTable<-read.table("SingaporeDietGenusCutOffFiltered.txt", quote="",header=T,sep="\t",row.names=1)

#This is a file containing sample level descriptions of the fungal communities at a genus resolution for samples from Bali
BaliFungusTable<-read.table("BaliFungusGenusCutOffFiltered.txt", quote="",header=T,sep="\t",row.names=1)

#This is a file containing a sample level description of the diet at a genus resolution for samples from Bali
BaliDietTable<-read.table("BaliDietGenusCutOffFiltered.txt", quote="",header=T,sep="\t",row.names=1)

#The following code conducts a partial Mantel test for the diet and mycobiome of Bali
#A partial Mantel test can use a third matrix to control for another variable
#Here we will use read count totals as our third matrix

#The following code calculates total reads for samples from Bali
BaliReadTotals<-colSums(rbind(BaliDietTable,BaliFungusTable))

#Here we set the seed to get consistent results
set.seed(6)
#The following code conducts a partial mantel test looking at correlations between the diet and mycobiome of samples from Bali while controlling for the effect of read count
#We are using Bray-Curtis distance and 100,000 permutations
ecodist::mantel(formula=(vegdist(t(BaliDietTable), method="bray", binary=FALSE)) ~(vegdist(t(BaliFungusTable), method="bray", binary=FALSE)) +(vegdist(BaliReadTotals, method="euclid")), nperm=100000)

#We focus on pval3, which is statistically significant
#We use this p-value because it is associated with a two-tailed test, which has a null hypothesis of r = 0

#The following code conducts a partial Mantel test for the diet and mycobiome of Singapore
#The following code calculates total reads for samples from Singapore
SingaporeReadTotals<-colSums(rbind(SingaporeDietTable,SingaporeFungusTable))

#Here we set the seed to get consistent results
set.seed(6)
#The following code conducts a partial mantel test looking at correlations between the diet and mycobiome of samples from Bali while controlling for the effect of read count
#We are using Bray-Curtis distance and 100,000 permutations
ecodist::mantel(formula=(vegdist(t(SingaporeDietTable), method="bray", binary=FALSE)) ~(vegdist(t(SingaporeFungusTable), method="bray", binary=FALSE)) +(vegdist(SingaporeReadTotals, method="euclid")), nperm=100000)

#This partial Mantel test is statistically significant