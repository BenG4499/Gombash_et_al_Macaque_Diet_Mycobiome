#Run this after SingaporeDietMycoMFA.R

#This is the cleanplot.pca function from Borcard et al 2018
'cleanplot.pca' <-
  function(res.pca,
           ax1 = 1,
           ax2 = 2,
           scaling = 1,
           plot.sites = TRUE,
           plot.spe = TRUE,
           label.sites = TRUE,
           label.spe = TRUE,
           cex.char1 = 0.7,
           pos.sites = 2,
           pos.spe = 4,
           mult.spe = 1,
           select.spe = NULL,
           mar.percent = 0.1,
           optimum = TRUE,
           move.origin = c(0, 0),
           silent = TRUE) {
    
    ### Internal functions
    'stretch' <-
      function(sites, mat, ax1, ax2, n, silent = silent) {
        # Compute stretching factor for the species arrows
        # First, compute the longest distance to centroid for the sites
        tmp1 <- rbind(c(0, 0), sites[, c(ax1, ax2)])
        D <- dist(tmp1)
        target <- max(D[1:n])
        # Then, compute the longest distance to centroid for the species arrows
        if (is.matrix(mat)) {
          p <- nrow(mat)   # Number of species to be drawn
          tmp2 <- rbind(c(0, 0), mat[, c(ax1, ax2)])
          D <- dist(tmp2)
          longest <- max(D[1:p])
        } else {
          tmp2 <- rbind(c(0, 0), mat[c(ax1, ax2)])
          longest <- dist(tmp2)
          # print(tmp2)
        }  # If a single row left in 'mat'
        #
        if (!silent)
          cat("target =",
              target,
              " longest =",
              longest,
              " fact =",
              target / longest,
              "\n")
        fact <- target / longest
      }
    
    'larger.plot' <-
      function(sit.sc,
               spe.sc,
               percent,
               move.origin,
               ax1,
               ax2) {
        # Internal function to expand plot limits (adapted from code by Pierre
        # Legendre)
        mat <- rbind(sit.sc, spe.sc)
        range.mat <- apply(mat, 2, range)
        rownames(range.mat) <- c("Min", "Max")
        z <- apply(range.mat, 2, function(x)
          x[2] - x[1])
        range.mat[1,] <- range.mat[1,] - z * percent
        range.mat[2,] <- range.mat[2,] + z * percent
        if (move.origin[1] != 0)
          range.mat[, ax1] <- range.mat[, ax1] - move.origin[1]
        if (move.origin[2] != 0)
          range.mat[, ax2] <- range.mat[, ax2] - move.origin[2]
        range.mat
      }
    
    #Here is the pcacircle function from Borcard et al 2018, it may be built into cleanplot.
    "pcacircle" <-
      function (pca, mult.spe, fact.spe, silent = silent) {
        # This function draws a circle of equilibrium contribution on a PCA plot
        # generated from the result file of a vegan rda() analysis.
        eigenv <- pca$CA$eig
        p <- length(eigenv)
        n <- nrow(pca$CA$u)
        tot <- sum(eigenv)
        radius <- (2 / p) ^ 0.5 * mult.spe * fact.spe
        symbols(
          0,
          0,
          circles = radius,
          inches = FALSE,
          add = TRUE,
          fg = 2
        )
        if (!silent) {
          cat(
            "\nSpecies arrows and the radius of the equilibrium circle are stretched ",
            "by a factor of",
            mult.spe * fact.spe
          )
          cat(
            "\nThe radius of the equilibrium circle is thus",
            (2 / p) ^ 0.5,
            "*",
            mult.spe,
            "*",
            fact.spe,
            "=",
            radius,
            "\n"
          )
        }
      }
    ### End internal functions
    
    if (!class(res.pca)[1] == "rda")
      stop("The input file is not a vegan output object of class 'rda'",
           call. = FALSE)
    if (!(is.null(res.pca$CCA)))
      stop(
        "The input file contains an RDA, not a PCA result. ",
        "Use function triplot.rda from the NEwR (2018) book to produce an RDA triplot."
      )
    if (scaling != 1 &
        scaling != 2)
      stop("Function only available for scaling 1 or 2", call. = FALSE)
    
    k <- length(res.pca$CA$eig)         # n. of PCA eigenvalues
    n.sp <- length(res.pca$colsum)      # n. of species
    ahead <- 0.05   # Length of arrow heads
    aangle <- 30    # Angle of arrow heads
    # 'vec' will contain the selection of species to be drawn
    if (is.null(select.spe)) {
      vec <- 1:n.sp
    } else {
      vec <- select.spe
    }
    
    # Scaling 1: the species scores have norms of 1
    # Scaling 1: the site scores are scaled to variances = can.eigenvalues
    # Scaling 2: the species scores have norms of sqrt(can.eigenvalues)
    # Scaling 2: the site scores are scaled to variances of 1
    
    # This version reconstructs and uses the original RDA output of L&L 2012,
    # Section 11.1.3
    
    Tot.var = res.pca$tot.chi         # Total variance in response data Y
    eig.val = res.pca$CA$eig          # Eigenvalues of Y-hat
    Lambda = diag(eig.val)            # Diagonal matrix of eigenvalues
    eig.val.rel = eig.val / Tot.var   # Relative eigenvalues of Y-hat
    Diag = diag(sqrt(eig.val.rel))    # Diagonal matrix of sqrt(relative
    # eigenvalues)
    U.sc1 = res.pca$CA$v              # Species scores, scaling=1
    U.sc2 = U.sc1 %*% Lambda ^ (0.5)  # Species scores, scaling=2
    n = nrow(res.pca$CA$u)            # Number of observations
    Z.sc2 = res.pca$CA$u * sqrt(n - 1)# Site scores, scaling=2
    Z.sc1 = Z.sc2 %*% Lambda ^ (0.5)  # Site scores, scaling=1
    
    if (is.null(select.spe)) {
      vec <- 1:n.sp
    } else {
      vec <- select.spe
    }
    
    if (scaling == 1) {
      sit.sc <- Z.sc1
      spe.sc <- U.sc1[vec,]
    } else {
      # For scaling=2
      sit.sc <- Z.sc2
      spe.sc <- U.sc2[vec,]
    }
    if (is.null(rownames(sit.sc)))
      rownames(sit.sc) <- paste("Site", 1:n, sep = "")
    if (is.null(rownames(spe.sc)))
      rownames(spe.sc) <- paste("Sp", 1:n.sp, sep = "")
    
    fact.spe <- 1
    if (optimum) {
      fact.spe <-
        stretch(sit.sc[, 1:k], spe.sc[, 1:k], ax1, ax2, n, silent = silent)
    }
    if (!silent)
      cat("fact.spe =", fact.spe, "\n\n")
    spe.sc <- spe.sc * fact.spe * mult.spe
    
    lim <-
      larger.plot(
        sit.sc[, 1:k],
        spe.sc[, 1:k],
        percent = mar.percent,
        move.origin = move.origin,
        ax1 = ax1,
        ax2 = ax2
      )
    if (!silent)
      print(lim)
    
    # Draw the main plot
    mat <- rbind(sit.sc[, 1:k], spe.sc[, 1:k])
    plot(
      mat[, c(ax1, ax2)],
      type = "n",
      main = paste("PCA biplot - Scaling", scaling),
      xlim = c(lim[1, ax1], lim[2, ax1]),
      ylim = c(lim[1, ax2], lim[2, ax2]),
      xlab = paste("PCA ", ax1),
      ylab = paste("PCA ", ax2),
      asp = 1
    )
    abline(h = 0, v = 0, col = "grey60")
    
    # Draw the site scores
    if (plot.sites) {
      points(sit.sc[, ax1], sit.sc[, ax2], pch = 20)
      if (label.sites)
        text(
          sit.sc[, ax1],
          sit.sc[, ax2],
          labels = rownames(sit.sc),
          col = "black",
          pos = pos.sites,
          cex = cex.char1
        )
    } else {
      if (label.sites)
        text(
          sit.sc[, ax1],
          sit.sc[, ax2],
          labels = rownames(sit.sc),
          col = "black",
          pos = NULL,
          cex = cex.char1
        )
    }
    
    # Draw the species scores
    if (plot.spe) {
      arrows(
        0,
        0,
        spe.sc[, ax1],
        spe.sc[, ax2],
        length = ahead,
        angle = aangle,
        col = "red"
      )
      if (label.spe)
        text(
          spe.sc[, ax1],
          spe.sc[, ax2],
          labels = rownames(spe.sc),
          col = "red",
          pos = pos.spe,
          cex = cex.char1
        )
    } else {
      if (label.spe)
        text(
          spe.sc[, ax1],
          spe.sc[, ax2],
          labels = rownames(spe.sc),
          col = "red",
          pos = NULL,
          cex = cex.char1
        )
    }
    
    # If scaling = 1 draw circle of equilibrium contribution
    if (scaling == 1) {
      pcacircle(
        res.pca,
        mult.spe = mult.spe,
        fact.spe = fact.spe,
        silent = silent
      )
    }
  }

DietFungiMFA$eig
#setwd("D:/01 Documents/School/ND/Macaque Diet/Publications/Diet Mycobiome/KeyTaxa")

#write.csv(data.frame(sort(DietFungiMFA$quanti.var$contrib[,1], decreasing = TRUE)),"KeyTaxa1S.csv")
#write.csv(data.frame(sort(DietFungiMFA$quanti.var$contrib[,2], decreasing = TRUE)),"KeyTaxa2S.csv")
#write.csv(data.frame(sort(DietFungiMFA$quanti.var$contrib[,3], decreasing = TRUE)),"KeyTaxa3S.csv")
#write.csv(data.frame(sort(DietFungiMFA$quanti.var$contrib[,4], decreasing = TRUE)),"KeyTaxa4S.csv")
#write.csv(data.frame(sort(DietFungiMFA$quanti.var$contrib[,5], decreasing = TRUE)),"KeyTaxa5S.csv")

fviz_contrib(DietFungiMFA,choice="quanti.var",axes=1,top=13,palette=c("#00FFFF","#A9A9A9","#FFA500")) +
  theme(text = element_text(size = 20))
fviz_contrib(DietFungiMFA,choice="quanti.var",axes=2,top=16,palette=c("#B8860B","#FF00FF","#FF0000","#00FFFF","#EE82EE","#A9A9A9","#FFA500")) +
  theme(text = element_text(size = 20))
fviz_contrib(DietFungiMFA,choice="quanti.var",axes=3,top=9,palette=c("#00FFFF","#EE82EE","#A9A9A9","#FFA500")) + theme(text = element_text(size = 20))
fviz_contrib(DietFungiMFA,choice="quanti.var",axes=4,top=26,palette=c("#FF00FF","#FF0000","#00FFFF","#EE82EE","#00CD00","#A9A9A9","#FFA500")) +
  theme(text = element_text(size = 20))
fviz_contrib(DietFungiMFA,choice="quanti.var",axes=5,top=22,palette=c("#B8860B","#FF00FF","#FF0000","#00FFFF","#EE82EE","#00CD00","#A9A9A9","#FFA500")) + theme(text = element_text(size = 20))


#Now I drop all of the 87 taxa from their tables
AnPathTable2 <- SingaporeFungusTable[which(SingaporeFungusKey[,2] == "Animal Pathogen" ),]
dim(AnPathTable2)
AnPathDrop<-c(which(rownames(AnPathTable2)=="Coccidioides"),
              which(rownames(AnPathTable2)=="Furia"))
AnPathTable2 <- AnPathTable2[-AnPathDrop,]
dim(AnPathTable2)
AnPathTable2<-t(AnPathTable2)
AnPathTable2<-decostand(AnPathTable2,method="hellinger")
AnPathTable2<-t(AnPathTable2)



PlPathTable2 <- SingaporeFungusTable[which(SingaporeFungusKey[,2] == "Plant Pathogen" ),]
dim(PlPathTable2)
PlPathDrop<-c(which(rownames(PlPathTable2)=="Asteroma"),
              which(rownames(PlPathTable2)=="Athelia"),
              which(rownames(PlPathTable2)=="Auriscalpium"),
              which(rownames(PlPathTable2)=="Botryosphaeria"),
              which(rownames(PlPathTable2)=="Caliciopsis"),
              which(rownames(PlPathTable2)=="Coprinus"),
              which(rownames(PlPathTable2)=="Eocronartium"),
              which(rownames(PlPathTable2)=="Limonomyces"),
              which(rownames(PlPathTable2)=="Macrophomina"),
              which(rownames(PlPathTable2)=="Passalora"),
              which(rownames(PlPathTable2)=="Pluteus"),
              which(rownames(PlPathTable2)=="Protomyces"),
              which(rownames(PlPathTable2)=="Schizonella"),
              which(rownames(PlPathTable2)=="Scleroconidioma"),
              which(rownames(PlPathTable2)=="Sympodiomycopsis"),
              which(rownames(PlPathTable2)=="Tilletiaria"))
PlPathTable2 <- PlPathTable2[-PlPathDrop,]
dim(PlPathTable2)
PlPathTable2<-t(PlPathTable2)
PlPathTable2<-decostand(PlPathTable2,method="hellinger")
PlPathTable2<-t(PlPathTable2)

LichTable2 <- SingaporeFungusTable[which(SingaporeFungusKey[,2] == "Lichenized" ),]
dim(LichTable2)
LichDrop<-c(which(rownames(LichTable2)=="Cyphellostereum"),
            which(rownames(LichTable2)=="Dictyocatenulata"),
            which(rownames(LichTable2)=="Digitothyrea"),
            which(rownames(LichTable2)=="Graphis"),
            which(rownames(LichTable2)=="Letrouitia"),
            which(rownames(LichTable2)=="Peccania"),
            which(rownames(LichTable2)=="Physconia"),
            which(rownames(LichTable2)=="Placopsis"),
            which(rownames(LichTable2)=="Rhizocarpon"),
            which(rownames(LichTable2)=="Sphaerulina"),
            which(rownames(LichTable2)=="Umbilicaria"))
LichTable2 <- LichTable2[-LichDrop,]
dim(LichTable2)
LichTable2<-t(LichTable2)
LichTable2<-decostand(LichTable2,method="hellinger")
LichTable2<-t(LichTable2)


PlASSTable2 <- SingaporeFungusTable[which(SingaporeFungusKey[,2] == "Plant Associated" ),]
dim(PlASSTable2)
PlAssDrop<-c(which(rownames(PlASSTable2)=="Ambispora"),
             which(rownames(PlASSTable2)=="Dentiscutata"),
             which(rownames(PlASSTable2)=="Gomphidius"),
             which(rownames(PlASSTable2)=="Gomphus"),
             which(rownames(PlASSTable2)=="Haplotrichum"),
             which(rownames(PlASSTable2)=="Hydnellum"),
             which(rownames(PlASSTable2)=="Hydnum"),
             which(rownames(PlASSTable2)=="Mortierella"),
             which(rownames(PlASSTable2)=="Phialocephala"),
             which(rownames(PlASSTable2)=="Rhizophagus"),
             which(rownames(PlASSTable2)=="Thelebolus"),
             which(rownames(PlASSTable2)=="Tomentella"),
             which(rownames(PlASSTable2)=="Troposporella"),
             which(rownames(PlASSTable2)=="Williopsis"))
PlASSTable2 <- PlASSTable2[-PlAssDrop,]
dim(PlASSTable2)
PlASSTable2<-t(PlASSTable2)
PlASSTable2<-decostand(PlASSTable2,method="hellinger")
PlASSTable2<-t(PlASSTable2)

MiscFungTable2 <- SingaporeFungusTable[which(SingaporeFungusKey[,2] == "Misc Fungus" ),]
dim(MiscFungTable2)
MiscFungDrop<-c(which(rownames(MiscFungTable2)=="Bandonia"),
                which(rownames(MiscFungTable2)=="Golubevia"),
                which(rownames(MiscFungTable2)=="Kuraishia"),
                which(rownames(MiscFungTable2)=="Marchandiomyces"),
                which(rownames(MiscFungTable2)=="Microsporomyces"),
                which(rownames(MiscFungTable2)=="Pterula"),
                which(rownames(MiscFungTable2)=="Saccharomycopsis"),
                which(rownames(MiscFungTable2)=="Teratoramularia"))
MiscFungTable2 <- MiscFungTable2[-MiscFungDrop,]
dim(MiscFungTable2)
MiscFungTable2<-t(MiscFungTable2)
MiscFungTable2<-decostand(MiscFungTable2,method="hellinger")
MiscFungTable2<-t(MiscFungTable2)

NonCropTable2 <- SingaporeDietTable[which(SingaporeGeneraCropList[1:which(rownames(SingaporeDietTable)=="Pinus"),1] == 0 & SingaporeGeneraMedList[1:which(rownames(SingaporeDietTable)=="Pinus"),4] == 0  ),]
dim(NonCropTable2)
NonCropDrop<-c(which(rownames(NonCropTable2)=="Connarus"),
                which(rownames(NonCropTable2)=="Eichhornia"),
                which(rownames(NonCropTable2)=="Ledebouria"))
NonCropTable2 <- NonCropTable2[-NonCropDrop,]
dim(NonCropTable2)
NonCropTable2<-t(NonCropTable2)
NonCropTable2<-decostand(NonCropTable2,method="hellinger")
NonCropTable2<-t(NonCropTable2)

CropTable2 <- SingaporeDietTable[which(SingaporeGeneraCropList[,1] == 1 ),]
dim(CropTable2)
CropDrop<-c(which(rownames(CropTable2)=="Asparagus"),
            which(rownames(CropTable2)=="Brassica"),
            which(rownames(CropTable2)=="Chenopodium"),
            which(rownames(CropTable2)=="Cucumis"),
            which(rownames(CropTable2)=="Flacourtia"),
            which(rownames(CropTable2)=="Manihot"),
            which(rownames(CropTable2)=="Morus"),
            which(rownames(CropTable2)=="Musa"),
            which(rownames(CropTable2)=="Persea"),
            which(rownames(CropTable2)=="Zea"),
            which(rownames(CropTable2)=="Zizania"))
CropTable2 <- CropTable2[-CropDrop,]
dim(CropTable2)
CropTable2<-t(CropTable2)
CropTable2<-decostand(CropTable2,method="hellinger")
CropTable2<-t(CropTable2)

MedTable2 <- SingaporeDietTable[which(SingaporeGeneraCropList[,1] == 0 & SingaporeGeneraMedList[,4] == 1  ),]
dim(MedTable2)
MedDrop<-c(which(rownames(MedTable2)=="Actinidia"),
           which(rownames(MedTable2)=="Averrhoa"),
           which(rownames(MedTable2)=="Bowiea"),
           which(rownames(MedTable2)=="Calophyllum"),
           which(rownames(MedTable2)=="Carallia"),
           which(rownames(MedTable2)=="Caryota"),
           which(rownames(MedTable2)=="Coptis"),
           which(rownames(MedTable2)=="Datisca"),
           which(rownames(MedTable2)=="Eucalyptus"),
           which(rownames(MedTable2)=="Gloriosa"),
           which(rownames(MedTable2)=="Humiria"),
           which(rownames(MedTable2)=="Larrea"),
           which(rownames(MedTable2)=="Licania"),
           which(rownames(MedTable2)=="Mazus"),
           which(rownames(MedTable2)=="Nepeta"),
           which(rownames(MedTable2)=="Osyris"),
           which(rownames(MedTable2)=="Posidonia"),
           which(rownames(MedTable2)=="Rhazya"),
           which(rownames(MedTable2)=="Ricinus"),
           which(rownames(MedTable2)=="Saruma"),
           which(rownames(MedTable2)=="Senna"))
MedTable2 <- MedTable2[-MedDrop,]
dim(MedTable2)
MedTable2<-t(MedTable2)
MedTable2<-decostand(MedTable2,method="hellinger")
MedTable2<-t(MedTable2)


#Here are the clean plot figures
#One to ID the key taxa, one to look nice

#Animal Pathogens
cleanplot.pca(rda(t(AnPathTable2)),scaling = 1, pos.spe=c(1,2,3),plot.spe=TRUE,label.sites=FALSE,label.spe=TRUE,select.spe=NULL,plot.sites=FALSE)
rownames(AnPathTable2)[-c(10,22,36,49)]<-""
cleanplot.pca(rda(t(AnPathTable2)),scaling = 1, pos.spe=c(1,2,3), plot.spe=TRUE,label.sites=FALSE,label.spe=TRUE,select.spe=NULL,plot.sites=FALSE)
#Plant Pathogens
cleanplot.pca(rda(t(PlPathTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
rownames(PlPathTable2)[-c(29,44,58)]<-""
cleanplot.pca(rda(t(PlPathTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
#Plant Associated
cleanplot.pca(rda(t(PlASSTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE,pos.spe=c(1,2,3))
rownames(PlASSTable2)[-c(8,17,18)]<-""
cleanplot.pca(rda(t(PlASSTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE,pos.spe=c(1,2,3))
#Saprotrophs
cleanplot.pca(rda(t(SapTable)),scaling = 1,plot.sites=FALSE,label.sites=FALSE, pos.spe=c(2,1))
rownames(SapTable)[-c(8,34,55,56,58,65)]<-""
cleanplot.pca(rda(t(SapTable)),scaling = 1,plot.sites=FALSE,label.sites=FALSE, pos.spe=c(1,2))
#Lichens
cleanplot.pca(rda(t(LichTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
rownames(LichTable2)[-c(9,16,17)]<-""
cleanplot.pca(rda(t(LichTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
#Miscellaneous fungi
cleanplot.pca(rda(t(MiscFungTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE, pos.spe=c(4,1))
rownames(MiscFungTable2)[-c(11,29,36,41,52)]<-""
cleanplot.pca(rda(t(MiscFungTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
#Crops
cleanplot.pca(rda(t(CropTable2)),scaling = 1,pos.spe=c(1,3),plot.sites=FALSE,label.sites=FALSE)
rownames(CropTable2)[-c(1,4,10,25,44)]<-""
cleanplot.pca(rda(t(CropTable2)),scaling = 1,pos.spe=c(1,3),plot.sites=FALSE,label.sites=FALSE)
#Anti-fungal plants
cleanplot.pca(rda(t(MedTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE, pos.spe=c(1,2))
rownames(MedTable2)[-c(3,13,26,39,44,45,62,63)]<-""
cleanplot.pca(rda(t(MedTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
#Miscellaneous plants
cleanplot.pca(rda(t(NonCropTable2)),scaling = 1,pos.spe=c(1,3),plot.sites=FALSE,label.sites=FALSE)
rownames(NonCropTable2)[-c(11,166,195,206,229)]<-""
cleanplot.pca(rda(t(NonCropTable2)),scaling = 1,pos.spe=c(4,1),plot.sites=FALSE,label.sites=FALSE)
#Invertebrates
cleanplot.pca(rda(t(ArthTable)),scaling = 1,plot.sites=FALSE,label.sites=FALSE,pos.spe=c(1,3))
rownames(ArthTable)[-c(34,36,98,107,132,141)]<-""
cleanplot.pca(rda(t(ArthTable)),scaling = 1,plot.sites=FALSE,label.sites=FALSE,pos.spe=c(2,3,4))
