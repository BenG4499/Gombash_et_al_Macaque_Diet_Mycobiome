#Run this after BaliDietMycoMFA.R

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

#write.csv(data.frame(sort(DietFungiMFA$quanti.var$contrib[,1], decreasing = TRUE)),"KeyTaxa1B.csv")
#write.csv(data.frame(sort(DietFungiMFA$quanti.var$contrib[,2], decreasing = TRUE)),"KeyTaxa2B.csv")
#write.csv(data.frame(sort(DietFungiMFA$quanti.var$contrib[,3], decreasing = TRUE)),"KeyTaxa3B.csv")
#write.csv(data.frame(sort(DietFungiMFA$quanti.var$contrib[,4], decreasing = TRUE)),"KeyTaxa4B.csv")
#write.csv(data.frame(sort(DietFungiMFA$quanti.var$contrib[,5], decreasing = TRUE)),"KeyTaxa5B.csv")

fviz_contrib(DietFungiMFA,choice="quanti.var",axes=1,top=24,palette=c("#B8860B","#FF00FF","#00FFFF","#EE82EE","#A9A9A9","#FFA500")) +
  theme(text = element_text(size = 20))
fviz_contrib(DietFungiMFA,choice="quanti.var",axes=2,top=23,palette=c("#B8860B","#FF00FF","#00FFFF","#EE82EE","#A9A9A9","#FFA500")) +
  theme(text = element_text(size = 20))
fviz_contrib(DietFungiMFA,choice="quanti.var",axes=3,top=21,palette=c("#B8860B","#FF00FF","#FF0000","#00FFFF","#EE82EE","#A9A9A9","#FFA500")) + theme(text = element_text(size = 20))
fviz_contrib(DietFungiMFA,choice="quanti.var",axes=4,top=19,palette=c("#B8860B","#FF00FF","#FF0000","#00FFFF","#EE82EE","#A9A9A9")) +
  theme(text = element_text(size = 20))
fviz_contrib(DietFungiMFA,choice="quanti.var",axes=5,top=29,palette=c("#B8860B","#FF00FF","#FF0000","#00FFFF","#EE82EE","#A9A9A9","#FFA500","#FFB6C1")) + theme(text = element_text(size = 20))


#Now I drop all of the 110 taxa from their tables
AnPathTable2 <- BaliFungusTable[which(BaliFungusKey[,2] == "Animal Pathogen" ),]
dim(AnPathTable2)
AnPathDrop<-c(which(rownames(AnPathTable2)=="Cladosporium"),
              which(rownames(AnPathTable2)=="Cochlonema"),
              which(rownames(AnPathTable2)=="Conidiobolus"),
              which(rownames(AnPathTable2)=="Coniosporium"),
              which(rownames(AnPathTable2)=="Entomophthora"),
              which(rownames(AnPathTable2)=="Fusarium"),
              which(rownames(AnPathTable2)=="Pyrenochaeta"),
              which(rownames(AnPathTable2)=="Schwanniomyces"),
              which(rownames(AnPathTable2)=="Scolecobasidium"),
              which(rownames(AnPathTable2)=="Sterigmatomyces"),
              which(rownames(AnPathTable2)=="Symbiotaphrina"),
              which(rownames(AnPathTable2)=="Tolypocladium"),
              which(rownames(AnPathTable2)=="Trigonopsis"),
              which(rownames(AnPathTable2)=="Verruconis"),
              which(rownames(AnPathTable2)=="Verticillium"))
AnPathTable2 <- AnPathTable2[-AnPathDrop,]
dim(AnPathTable2)
AnPathTable2<-t(AnPathTable2)
AnPathTable2<-decostand(AnPathTable2,method="hellinger")
AnPathTable2<-t(AnPathTable2)



PlPathTable2 <- BaliFungusTable[which(BaliFungusKey[,2] == "Plant Pathogen" ),]
dim(PlPathTable2)
PlPathDrop<-c(which(rownames(PlPathTable2)=="Acremonium"),
              which(rownames(PlPathTable2)=="Acrodontium"),
              which(rownames(PlPathTable2)=="Colletotrichum"),
              which(rownames(PlPathTable2)=="Fomitiporia"),
              which(rownames(PlPathTable2)=="Itersonilia"),
              which(rownames(PlPathTable2)=="Lasiodiplodia"),
              which(rownames(PlPathTable2)=="Leucocintractia"),
              which(rownames(PlPathTable2)=="Marssonina"),
              which(rownames(PlPathTable2)=="Olpidium"),
              which(rownames(PlPathTable2)=="Plectosphaerella"),
              which(rownames(PlPathTable2)=="Pluteus"),
              which(rownames(PlPathTable2)=="Rickenella"),
              which(rownames(PlPathTable2)=="Sclerotium"),
              which(rownames(PlPathTable2)=="Truncatella"),
              which(rownames(PlPathTable2)=="Ustilentyloma"))
PlPathTable2 <- PlPathTable2[-PlPathDrop,]
dim(PlPathTable2)
PlPathTable2<-t(PlPathTable2)
PlPathTable2<-decostand(PlPathTable2,method="hellinger")
PlPathTable2<-t(PlPathTable2)

LichTable2 <- BaliFungusTable[which(BaliFungusKey[,2] == "Lichenized" ),]
dim(LichTable2)
LichDrop<-c(which(rownames(LichTable2)=="Anamylopsora"),
            which(rownames(LichTable2)=="Anzia"),
            which(rownames(LichTable2)=="Combea"),
            which(rownames(LichTable2)=="Hypogymnia"),
            which(rownames(LichTable2)=="Lecidea"),
            which(rownames(LichTable2)=="Lopezaria"),
            which(rownames(LichTable2)=="Ochrolechia"),
            which(rownames(LichTable2)=="Pannoparmelia"),
            which(rownames(LichTable2)=="Placopsis"),
            which(rownames(LichTable2)=="Pyrgillus"),
            which(rownames(LichTable2)=="Sphaerulina"),
            which(rownames(LichTable2)=="Umbilicaria"))
LichTable2 <- LichTable2[-LichDrop,]
dim(LichTable2)
LichTable2<-t(LichTable2)
LichTable2<-decostand(LichTable2,method="hellinger")
LichTable2<-t(LichTable2)


PlASSTable2 <- BaliFungusTable[which(BaliFungusKey[,2] == "Plant Associated" ),]
dim(PlASSTable2)
PlAssDrop<-c(which(rownames(PlASSTable2)=="Archaeospora"),
             which(rownames(PlASSTable2)=="Cenococcum"),
             which(rownames(PlASSTable2)=="Dioszegia"),
             which(rownames(PlASSTable2)=="Diversispora"),
             which(rownames(PlASSTable2)=="Elaphomyces"),
             which(rownames(PlASSTable2)=="Gigaspora"),
             which(rownames(PlASSTable2)=="Haplotrichum"),
             which(rownames(PlASSTable2)=="Hygrophorus"),
             which(rownames(PlASSTable2)=="Inocybe"),
             which(rownames(PlASSTable2)=="Kurtzmanomyces"),
             which(rownames(PlASSTable2)=="Papiliotrema"),
             which(rownames(PlASSTable2)=="Paxillus"),
             which(rownames(PlASSTable2)=="Rhizophagus"),
             which(rownames(PlASSTable2)=="Rhizopogon"),
             which(rownames(PlASSTable2)=="Strobiloscypha"),
             which(rownames(PlASSTable2)=="Tricholoma"),
             which(rownames(PlASSTable2)=="Truncocolumella"),
             which(rownames(PlASSTable2)=="Vishniacozyma"),
             which(rownames(PlASSTable2)=="Williopsis"))
PlASSTable2 <- PlASSTable2[-PlAssDrop,]
dim(PlASSTable2)
PlASSTable2<-t(PlASSTable2)
PlASSTable2<-decostand(PlASSTable2,method="hellinger")
PlASSTable2<-t(PlASSTable2)

SapTable2 <- BaliFungusTable[which(BaliFungusKey[,2] == "Saprotroph" ),]
dim(SapTable2)
SapDrop<-c(which(rownames(SapTable2)=="Buckleyzyma"))
SapTable2 <- SapTable2[-SapDrop,]
dim(SapTable2)
SapTable2<-t(SapTable2)
SapTable2<-decostand(SapTable2,method="hellinger")
SapTable2<-t(SapTable2)

MiscFungTable2 <- BaliFungusTable[which(BaliFungusKey[,2] == "Misc Fungus" ),]
dim(MiscFungTable2)
MiscFungDrop<-c(which(rownames(MiscFungTable2)=="Classicula"),
                which(rownames(MiscFungTable2)=="Erythrobasidium"),
                which(rownames(MiscFungTable2)=="Friedmanniomyces"),
                which(rownames(MiscFungTable2)=="Golubevia"),
                which(rownames(MiscFungTable2)=="Kuraishia"),
                which(rownames(MiscFungTable2)=="Nais"),
                which(rownames(MiscFungTable2)=="Orbimyces"),
                which(rownames(MiscFungTable2)=="Piskurozyma"),
                which(rownames(MiscFungTable2)=="Pseudobensingtonia"),
                which(rownames(MiscFungTable2)=="Pseudohyphozyma"),
                which(rownames(MiscFungTable2)=="Pterula"),
                which(rownames(MiscFungTable2)=="Richoniella"),
                which(rownames(MiscFungTable2)=="Slooffia"),
                which(rownames(MiscFungTable2)=="Striaticonidium"),
                which(rownames(MiscFungTable2)=="Sugiyamaella"),
                which(rownames(MiscFungTable2)=="Unknown"))
MiscFungTable2 <- MiscFungTable2[-MiscFungDrop,]
dim(MiscFungTable2)
MiscFungTable2<-t(MiscFungTable2)
MiscFungTable2<-decostand(MiscFungTable2,method="hellinger")
MiscFungTable2<-t(MiscFungTable2)

CropTable2 <- BaliDietTable[which(BaliGeneraCropList[,1] == 1 ),]
dim(CropTable2)
CropDrop<-c(which(rownames(CropTable2)=="Aphandra"),
            which(rownames(CropTable2)=="Curcuma"),
            which(rownames(CropTable2)=="Elaeis"),
            which(rownames(CropTable2)=="Flacourtia"),
            which(rownames(CropTable2)=="Fragaria"),
            which(rownames(CropTable2)=="Helianthus"),
            which(rownames(CropTable2)=="Hordeum"),
            which(rownames(CropTable2)=="Ipomoea"),
            which(rownames(CropTable2)=="Persea"),
            which(rownames(CropTable2)=="Polianthes"),
            which(rownames(CropTable2)=="Rubus"),
            which(rownames(CropTable2)=="Saccharum"),
            which(rownames(CropTable2)=="Sorghum"),
            which(rownames(CropTable2)=="Triticum"),
            which(rownames(CropTable2)=="Zea"))
CropTable2 <- CropTable2[-CropDrop,]
dim(CropTable2)
CropTable2<-t(CropTable2)
CropTable2<-decostand(CropTable2,method="hellinger")
CropTable2<-t(CropTable2)

MedTable2 <- BaliDietTable[which(BaliGeneraCropList[,1] == 0 & BaliGeneraMedList[,4] == 1  ),]
dim(MedTable2)
MedDrop<-c(which(rownames(MedTable2)=="Arabidopsis"),
           which(rownames(MedTable2)=="Bowiea"),
           which(rownames(MedTable2)=="Callirhoe"),
           which(rownames(MedTable2)=="Calotropis"),
           which(rownames(MedTable2)=="Cordyline"),
           which(rownames(MedTable2)=="Dahlia"),
           which(rownames(MedTable2)=="Ficus"),
           which(rownames(MedTable2)=="Geranium"),
           which(rownames(MedTable2)=="Hedychium"),
           which(rownames(MedTable2)=="Larrea"),
           which(rownames(MedTable2)=="Monarda"),
           which(rownames(MedTable2)=="Neolamarckia"),
           which(rownames(MedTable2)=="Nothapodytes"),
           which(rownames(MedTable2)=="Oxalis"),
           which(rownames(MedTable2)=="Picea"),
           which(rownames(MedTable2)=="Rhazya"),
           which(rownames(MedTable2)=="Rosa"))
MedTable2 <- MedTable2[-MedDrop,]
dim(MedTable2)
MedTable2<-t(MedTable2)
MedTable2<-decostand(MedTable2,method="hellinger")
MedTable2<-t(MedTable2)


#Here are the clean plot figures
#One to ID the key taxa, one to look nice

#Animal Pathogens
cleanplot.pca(rda(t(AnPathTable2)),scaling = 1, plot.spe=TRUE,label.sites=FALSE,label.spe=TRUE,select.spe=NULL,plot.sites=FALSE)
rownames(AnPathTable2)[-c(10,24,33,57)]<-""
cleanplot.pca(rda(t(AnPathTable2)),scaling = 1, plot.spe=TRUE,label.sites=FALSE,label.spe=TRUE,select.spe=NULL,plot.sites=FALSE)
#Plant Pathogens
cleanplot.pca(rda(t(PlPathTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
rownames(PlPathTable2)[-c(4,29,43,78)]<-""
cleanplot.pca(rda(t(PlPathTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
#Plant Associated
cleanplot.pca(rda(t(PlASSTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
rownames(PlASSTable2)[-c(16,26)]<-""
cleanplot.pca(rda(t(PlASSTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
#Saprotrophs
cleanplot.pca(rda(t(SapTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
rownames(SapTable2)[-c(13,30,64,67)]<-""
cleanplot.pca(rda(t(SapTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
#Lichens
cleanplot.pca(rda(t(LichTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
rownames(LichTable2)[-c(13,19)]<-""
cleanplot.pca(rda(t(LichTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
#Miscellaneous fungi
cleanplot.pca(rda(t(MiscFungTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
rownames(MiscFungTable2)[-c(29,37,50,51,66)]<-""
cleanplot.pca(rda(t(MiscFungTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
#Crops
cleanplot.pca(rda(t(CropTable2)),scaling = 1,pos.spe=c(1,3),plot.sites=FALSE,label.sites=FALSE)
rownames(CropTable2)[-c(9,24,41,42)]<-""
cleanplot.pca(rda(t(CropTable2)),scaling = 1,pos.spe=c(1,3),plot.sites=FALSE,label.sites=FALSE)
#Anti-fungal plants
cleanplot.pca(rda(t(MedTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
rownames(MedTable2)[-c(40,48,63,65,66)]<-""
cleanplot.pca(rda(t(MedTable2)),scaling = 1,plot.sites=FALSE,label.sites=FALSE)
#Miscellaneous plants
cleanplot.pca(rda(t(NonCropTable)),scaling = 1,pos.spe=c(4,1),plot.sites=FALSE,label.sites=FALSE)
rownames(NonCropTable)[-c(59,110,143,144,187,194,299,313)]<-""
cleanplot.pca(rda(t(NonCropTable)),scaling = 1,pos.spe=c(4,1),plot.sites=FALSE,label.sites=FALSE)
#Invertebrates
cleanplot.pca(rda(t(ArthTable)),scaling = 1,plot.sites=FALSE,label.sites=FALSE,pos.spe=c(3,4))
rownames(ArthTable)[-c(60,115,136,158)]<-""
cleanplot.pca(rda(t(ArthTable)),scaling = 1,plot.sites=FALSE,label.sites=FALSE,pos.spe=c(3,4))
