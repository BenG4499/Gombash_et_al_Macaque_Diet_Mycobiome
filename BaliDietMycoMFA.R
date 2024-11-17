#Clear the working environment
rm(list=ls())
#Load the required packages
library(FactoMineR)#Contains code for running a Multiple Factor Analysis (MFA)
library(factoextra)#Contains functions to visualize results
library(vegan)#Contains decostand function for Hellinger transformation

#Setting the working directory, this will vary by machine
setwd("D:/01 Documents/School/ND/Macaque Diet/Publications/Diet Mycobiome/Data")

#Loading in the various data files
#This is a file containing a list of dietary genera found in samples from Bali and whether they are crop plants or not
BaliGeneraCropList<-read.table("BaliGeneraCropList.txt", quote="",header=T,sep="\t",row.names=1)

#This is a file containing a list of dietary genera found in samples from Bali and whether they have medicinal properties
BaliGeneraMedList<-read.table("BaliGeneraMedList.txt", quote="",header=T,sep="\t",row.names=1)

#This is a file containing a list of fungal genera found in samples from Bali and their assigned fungal category
#Fungal categories were assigned based on FunGUILD v 1.1 and a literature search
#Genera that could be placed into multiple categories were assigned to the category with the highest priority (see manuscript)
BaliFungusKey<-data.frame(read.table("BaliFungusKey.txt", quote="",header=T,sep="\t",row.names=1,check.names=FALSE))

#This is a file containing sample level descriptions of the fungal communities at a genus resolution for samples from Bali
BaliFungusTable<-read.table("BaliFungusGenusCutOffFiltered.txt", quote="",header=T,sep="\t",row.names=1)

#This is a file containing a sample level description of the diet at a genus resolution for samples from Bali
BaliDietTable<-read.table("BaliDietGenusCutOffFiltered.txt", quote="",header=T,sep="\t",row.names=1)


#Multiple factor analysis starts with several distinct tables, each table represents a group
#A Principal Components Analysis (PCA) is run on each table
#The first eigenvalue for each PCA is used to standardize its table
#These standardized tables combine and a PCA is run on this table

#This code uses the Fungus Key to select all Animal Pathogens and put them into a table, which is Hellinger transformed
AnPathTable <- BaliFungusTable[which(BaliFungusKey[,2] == "Animal Pathogen"),]
dim(AnPathTable)
AnPathTableReads<-AnPathTable
AnPathTable<-t(AnPathTable)
AnPathTable<-decostand(AnPathTable,method="hellinger")
AnPathTable<-t(AnPathTable)

#This code uses the Fungus Key to select all Plant Pathogens and put them into a table, which is Hellinger transformed
PlPathTable <- BaliFungusTable[which(BaliFungusKey[,2] == "Plant Pathogen" ),]
dim(PlPathTable)
PlPathTableReads<-PlPathTable
PlPathTable<-t(PlPathTable)
PlPathTable<-decostand(PlPathTable,method="hellinger")
PlPathTable<-t(PlPathTable)

#This code uses the Fungus Key to select all Lichenized Fungi and put them into a table, which is Hellinger transformed
LichTable <- BaliFungusTable[which(BaliFungusKey[,2] == "Lichenized" ),]
dim(LichTable)
LichTableReads<-LichTable
LichTable<-t(LichTable)
LichTable<-decostand(LichTable,method="hellinger")
LichTable<-t(LichTable)

#This code uses the Fungus Key to select all Plant Associated Fungi and put them into a table, which is Hellinger transformed
PlASSTable <- BaliFungusTable[which(BaliFungusKey[,2] == "Plant Associated" ),]
dim(PlASSTable)
PlASSTableReads<-PlASSTable
PlASSTable<-t(PlASSTable)
PlASSTable<-decostand(PlASSTable,method="hellinger")
PlASSTable<-t(PlASSTable)

#This code uses the Fungus Key to select all Saprotrophic Fungi and put them into a table, which is Hellinger transformed
SapTable <- BaliFungusTable[which(BaliFungusKey[,2] == "Saprotroph" ),]
dim(SapTable)
SapTableReads<-SapTable
SapTable<-t(SapTable)
SapTable<-decostand(SapTable,method="hellinger")
SapTable<-t(SapTable)

#This code uses the Fungus Key to select all Miscellaneous Fungi and put them into a table, which is Hellinger transformed
MiscFungTable <- BaliFungusTable[which(BaliFungusKey[,2] == "Misc Fungus" ),]
dim(MiscFungTable)
MiscFungTableReads<-MiscFungTable
MiscFungTable<-t(MiscFungTable)
MiscFungTable<-decostand(MiscFungTable,method="hellinger")
MiscFungTable<-t(MiscFungTable)

#This code uses the Crop Genera List to select all Crop Plants and put them into a table, which is Hellinger transformed
#Note that the Crop group has higher priority than the Anti-Fungal Medicinal Plant group
CropTable <- BaliDietTable[which(BaliGeneraCropList[,1] == 1 ),]
dim(CropTable)
CropTableReads<-CropTable
CropTable<-t(CropTable)
CropTable<-decostand(CropTable,method="hellinger")
CropTable<-t(CropTable)

#This code selects all Invertebrates and puts them into a table, which is Hellinger transformed
#Note, "Pinus" is the last plant genus in our data set, which is likely not true for other data
MedTable <- BaliDietTable[which(BaliGeneraCropList[,1] == 0 & BaliGeneraMedList[,4] == 1  ),]
dim(MedTable)
MedTableReads<-MedTable
MedTable<-t(MedTable)
MedTable<-decostand(MedTable,method="hellinger")
MedTable<-t(MedTable)

#This code selects all Invertebrates and puts them into a table, which is Hellinger transformed
#Note, "Pinus" is the last plant genus in our data set, which is likely not true for other data
ArthTable <- BaliDietTable[(which(rownames(BaliDietTable)=="Pinus")+1):dim(BaliDietTable)[1],]
dim(ArthTable)
ArthTableReads<-ArthTable
ArthTable<-t(ArthTable)
ArthTable<-decostand(ArthTable,method="hellinger")
ArthTable<-t(ArthTable)

#This code uses the Medicinal Genera List and Crop Genera List to select all Miscellaneous Plants, which are not crops and do not have Anti-Fungal properties, and put them into a table, which is Hellinger transformed
NonCropTable <- BaliDietTable[which(BaliGeneraCropList[1:which(rownames(BaliDietTable)=="Pinus"),1] == 0 & BaliGeneraMedList[1:which(rownames(BaliDietTable)=="Pinus"),4] == 0  ),]
dim(NonCropTable)
NonCropTableReads<-NonCropTable
NonCropTable<-t(NonCropTable)
NonCropTable<-decostand(NonCropTable,method="hellinger")
NonCropTable<-t(NonCropTable)

#Here we combine all of the previously created tables to prepare for the Multiple Factor Analysis
ComboTable<-rbind(CropTable,NonCropTable,ArthTable,MedTable,AnPathTable,PlASSTable,PlPathTable,LichTable,MiscFungTable,SapTable)
dim(ComboTable)
ComboTable<-t(ComboTable)
dim(ComboTable)

#This code runs the Multiple Factor Analysis
DietFungiMFA<-MFA(base=ComboTable,group = 
                    c(dim(CropTable)[1],dim(NonCropTable)[1],dim(ArthTable)[1],dim(MedTable)[1],dim(AnPathTable)[1],dim(PlASSTable)[1],dim(PlPathTable)[1],dim(LichTable)[1],dim(MiscFungTable)[1],dim(SapTable)[1]),
                  name.group=c("Crops","MiscPlant","Inverts","AntiFung","AniPath","PlantAssoc","PlantPath","Lichens","MiscFung","Sapro"),graph=TRUE)

#This code creates a color coded plot of relationships between the groups that were included in the MFA along dimensions 1 and 2
par(mfrow=c(1,1),mar=c(5.1, 5.1, 4.1, 2.1))
plot(y=DietFungiMFA$group$coord[,2],x=DietFungiMFA$group$coord[,1],pch = 19, ylab= "Dimension 2 (5.26%)",xlab="Dimension 1 (6.52%)", main = "Groups Representation",cex.lab =2, cex=2,ylim=c(0,1),xlim=c(0,1),col=c("#FF0000","#00CD00","#0000FF","#FF00FF","#B8860B","#A9A9A9","#FFA500","#00FFFF","#EE82EE","#FFB6C1"))
text(y=DietFungiMFA$group$coord[,2],x=DietFungiMFA$group$coord[,1],labels = c("Crops","MiscPlant","Inverts","AntiFung","AniPath","PlantAssoc","PlantPath","Lichens","MiscFung","Sapro"), pos=c(1,1,1,1,3,3,1,3,1,4))

#This code creates a color coded bar chart showing how much variation each group contributes to the first dimension of the MFA
fviz_contrib(DietFungiMFA,"group",axes=1, fill = "name", color = "name") +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values=c("#FFB6C1","#A9A9A9","#EE82EE","#00CD00","#B8860B","#FFA500","#0000FF","#00FFFF","#FF00FF","#FF0000"), name="Groups") +
  scale_color_manual(values=c("#FFB6C1","#A9A9A9","#EE82EE","#00CD00","#B8860B","#FFA500","#0000FF","#00FFFF","#FF00FF","#FF0000"), name="Groups")

#This code creates a color coded bar chart showing how much variation each group contributes to the second dimension of the MFA
fviz_contrib(DietFungiMFA,"group",axes=2, fill = "name", color = "name") +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values=c("#B8860B","#FFA500","#A9A9A9","#FFB6C1","#EE82EE","#00CD00","#00FFFF","#0000FF","#FF00FF","#FF0000"), name="Groups") +
  scale_color_manual(values=c("#B8860B","#FFA500","#A9A9A9","#FFB6C1","#EE82EE","#00CD00","#00FFFF","#0000FF","#FF00FF","#FF0000"), name="Groups")

#This code creates a color coded bar chart showing how much variation the top twenty-eight factors/genera contribute to the first dimension of the MFA
fviz_contrib(DietFungiMFA,choice="quanti.var",axes=1,top=28,palette=c("#B8860B","#00FFFF","#EE82EE","#A9A9A9","#FFA500")) +
  theme(text = element_text(size = 20))

#This code creates a color coded bar chart showing how much variation the top eighteen factors/genera contribute to the second dimension of the MFA
fviz_contrib(DietFungiMFA,choice="quanti.var",axes=2,top=18,palette=c("#B8860B","#FF00FF","#00FFFF","#EE82EE","#A9A9A9","#FFA500"))
