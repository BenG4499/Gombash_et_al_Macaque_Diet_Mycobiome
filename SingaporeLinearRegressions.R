#Clear the working environment
rm(list=ls())

#Load the required packages
library(car)
library(glmnet)

#Setting the working directory, this will vary by machine
setwd("D:/01 Documents/School/ND/Macaque Diet/Publications/Diet Mycobiome/Data")

#Loading in the linear regression data file
#This is a file containing a list of list of samples and the richness values associated with those samples.
#The richness values of fungal groups from the Multiple Factor Analysis are mutually exclusive, which is why they sum to the overall mycobiome richness
#The richness values of diet groups from the Multiple Factor Analysis are not mutually exclusive due to the overlap between crops and anti-fungal plants
#Invertebrate richness is not included to avoid aliasing issues and because interactions between fungi and plant diet items are of greater interest to the authors
SingaporeData<-read.table("SingaporeLinearRegressionData.txt", quote="",header=T,sep="\t",row.names=1)

#The following code runs a linear regression looking for a relationship between mycobiome richness and diet richness
MycoVsDietModel<-lm(MycoRich~DietRich, data = SingaporeData)
summary(MycoVsDietModel)
#This is statistically significant

#The following code creates a plot of the previously created model
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
plot(y=SingaporeData$MycoRich,x=SingaporeData$DietRich,pch = 19, ylab= "Fungal Genera Richness",xlab="Diet Genera Richness",main="Singapore", cex.lab =2, cex=2)
abline(MycoVsDietModel,col="red", lwd=2)

#The remaining code looks at doing a multivariate multiple linear regression. This involves more than one response variable (multivariate) 
#and more than one predictor variable (multiple)
#This will address questions about whether specific plant diet groups have additional impacts on fungal groups beyond the overall diet richness

#This code looks at the variance inflation factors to see if there is a problematic amount of multicollinearity between the diet predictors
#when we are attempting to predict the richness values of the fungal groups
vif(lm(PlantAssocRich + PlantPathRich + AniPathRich + SaproRich + LichenRich + MiscFungRich ~ CropRich + PlantRich + AntiFungRich + DietRich, data = SingaporeData))
#We have a problematic amount of multicollinearity here, the VIF exceed "5 or 10" (James et al., 2014)

#To address the multicollinearity, we are going to employ LASSO to remove some of the predictors
#LASSO depends on a variable called lambda to eliminate predictors
#The following code is one process for selecting a lambda value

#The following code sets a seed for consistent results
set.seed(6)
#The following code runs a cross validation method to select a lambda value
#The "mgaussian" family is specific for multivariate multiple linear regression
cvmodel<-cv.glmnet(data.matrix(SingaporeData[,8:11]),y = data.matrix(SingaporeData[,1:6]), alpha = 1, family = "mgaussian")
cvmodel$lambda.min
cvmodel$lambda.1se
plot(cvmodel)
#The values that we see are ~0.87 and ~6.1


#This code tests multiple potential lambda values between 0 and 10, with the two values that were previously identified inserted
bestmodel <- glmnet(x = data.matrix(SingaporeData[,8:11]),y = data.matrix(SingaporeData[,1:6]),
                    alpha = 1, family = "mgaussian", lambda = c(10,9,8,7,6,cvmodel$lambda.1se,4,3,2,cvmodel$lambda.min,0))
#The output that we are interested in is which parameters are dropped at various lambda values
coef(bestmodel,s = c(10,9,8,6,cvmodel$lambda.1se,4,3,2,cvmodel$lambda.min,0))
#A lambda value of 0 does not drop any parameters
#A lambda value of ~0.87 drops crop richness, plant richness, and anti-fungal plant richness
#We don't have any values that preserve a predictor other than overall diet richness

