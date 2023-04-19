###############################
####
####      Analyse du comportement_Exp daphnies
####      Auteur-trice : Christophe Reis & Floriane Tisserand
####      Date : 04/04/2023
####
####
####      source : 
####      trajr : McLean DJ, Skowron Volponi MA. trajr: An R package for 
####              characterisation of animal trajectories. Ethology. 2018;00:19. 
####              https://doi.org/10.1111/eth.12739.
####      TrackDem : https://github.com/marjoleinbruijning/trackdem
####      av : https://CRAN.R-project.org/package=av
####
###############################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Library and WD ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(trackdem)
library(av)
library(trajr)

# defini le WD
setwd("path")

#~~~~~~~~~~~~~~~~~~~~ selection de l'experience ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ne runner que la ligne correspondantea l'experience
#exp = "low_conc"
exp = "high_conc_24h"
#exp = "high_conc_48h"

#~~~~~~~~~~~~~~~~~~ informations vidéos et DL vidéos ~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(exp == "low_conc"){ # videos low conc.

  conc <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10") # concentrations
  individus <- c(12,14,7,7,6,9,8,9,13,8) # nbre d'individus
  thres <- c(-0.12,-0.18,-0.16,-0.18,-0.21,-0.21,-0.115,-0.14,-0.18,-0.27) # seuil
  # telecharge les videos low conc.
  list_mov <- list.files(path = "videos/brute", pattern = (".mp4"), full.names = TRUE)
  
}

if(exp == "high_conc_24h"){ # videos high conc. 24h
  
  conc <- c("c0","c1","c2","c3","c4") # concentrations
  individus <- c(8,9,9,10,7) # nbre d'individus
  thres <- c(-0.115,-0.08,-0.115,-0.3,-0.3) # seuil
  # telecharge les videos high conc. 24h
  list_mov <- list.files(path = "videos/brute/24H", pattern = (".mp4"), full.names = TRUE)
  
}

if(exp == "high_conc_48h"){ # videos high conc. 48h
  
  conc <- c("c0","c1","c2","c3") # concentrations
  individus <- c(5,9,10,10) # nbre d'individus
  thres <- c(-0.115,-0.115,-0.07,-0.115) # seuil
  # telecharge les videos high conc. 48h 48H
  list_mov <- list.files(path = "videos/brute/48H", pattern = (".mp4"), full.names = TRUE)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ variables generales ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# nombre de FPS
FramePS <- 15

# creation tableau all_data pour sauvegarder les donnees
data_all = matrix(ncol = 12)
colnames(data_all) <- c("Concentration",
                        "distance [mm]",
                        "Vitesse moyenne [mm/s]",
                        "SD vitesse [mm/s]",
                        "accélération [mm/s^2]",
                        "Temps d'arret [s]", 
                        "straightness (D/L) [-]",
                        "alternative straigthness (r) [-]",
                        "sinuosite [-]",
                        "Emax [-]",
                        "DC [-]",
                        "SDDC [-]")

#~~~~~~~~~~~~~~~~~~~~~~~~~ remise à 0 du tableau des données ~~~~~~~~~~~~~~~~~~

file.remove(list.files(path = "resultats/imagesChemin", 
                       pattern = (".pdf"), 
                       full.names = TRUE))
file.remove(list.files(path = "resultats/videoChemin", 
                       pattern = (".mp4"), 
                       full.names = TRUE))
file.remove(list.files(path = "resultats/", 
                       pattern = (".csv"), 
                       full.names = TRUE))

#~~~~~~~~~~~~~~~~~~ Convert video to images sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (o in 1:length(conc)){ # commencement de la boucle

  # verification du depart d'une nouvelle boucle
  print(paste("On commence la concentration num :",o))
  
  # supprime les anciennes images du repertoires
  file.remove(list.files(path = "videos/images", 
                         pattern = (".png"), 
                         full.names = TRUE))

  # Coupe la video en images
  av_video_images(list_mov[o], 
                  destdir = "videos/images", 
                  format = "png", fps = FramePS)

  #~~~~~~~~~~~~~~~~~~~~~~ telechargement des images ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  compteur <- length(list.files("videos/images"))-1 # compte le nombre d'images
  N.images <- c(1:compteur) # definie un vct avec le nombre d'image
  allFullImages <- loadImages (dirPictures = "videos/images",
                               nImages = N.images)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~Suppression du bruit de fond~~~~~~~~~~~~~~~~~~~~~~

  #1) Creation de l'image bruit de fond
  stillBack <- createBackground(allFullImages,method="mean")

  #2) Suppression du fond 
  allImages <- subtractBackground(stillBack)

  #~~~~~~~~~~~~~~~~~~~ Identification des particules mobiles ~~~~~~~~~~~~~~~~~~~
  
  partIden <- identifyParticles(allImages,
                                threshold = thres[o],
                                pixelRange = c(5,500)) # reporter la valeur seuile
  
  print(summary(partIden)) # montre le resultat final pour controle
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~ Enregistrement des mouvements ~~~~~~~~~~~~~~~~~~~~~

  records <- trackParticles(partIden,L=200,R=3) # traque et lie les points. 

  #~~~~~~~~~~~~~~~~~~~ créer l'image de suivit et la video ~~~~~~~~~~~~~~~~~~~~~
  
  setwd("path")
  pdf(file = paste("imageSuivis_C",as.character(o),".pdf",sep=""))
  par(mfrow=c(1,1))
  plot(records,type="trajectories",incThres=1) # affichage du resultat
  dev.off()
  setwd("path")

  plot(records,type="animation",incThres=1,
       path="resultats/videoChemin/",
       libavpath = "libav/usr/bin/avconv.exe") # creation video de suivis
  file.rename(from="resultats/videoChemin/animation.mp4", 
              to = paste("resultats/videoChemin/",as.character(o),".mp4",sep="")) # renome le fichier
  
  #~~~~~~~~~~~~~~~~~~~~~~~ Enregistrement des coordonnées ~~~~~~~~~~~~~~~~~~~~~~
  
  coordonnes <- matrix(ncol = (individus[o]*2), nrow = (length(records$trackRecord[1,,1]))) # creation d'une matrice coordonnees
  
  wu <- 1
  wuu <- 1
  
  for(i in 1 : (individus[o]*2)){
    
    if(wuu==3){wu <- wu + 1}
    if(wuu==3){wuu <- 1}
    
    coordonnes[,i] <- records$trackRecord[wu,,wuu]
    
    wuu <- wuu + 1
    
  }
  
  if(o==1){data_coords = matrix(nrow = compteur)} # créer le tableau pour les coordonnees
  data_coords <- cbind(data_coords,coordonnes) # rajoute les donnees aux tableaux
  
  #~~~~~~~~~~~~~~~~~~~~~~~ Enregistrement des donnees ~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  resultat <- summary(records, incThres = 1)
  list_data <- Map(as.data.frame, resultat) # extrait les donnees
  
  #~~~~~~~ donnees utile pour les traits comportemantales ~~~~~~~~~~~~~~~~~~~~~~
  
  nbre <- summary(records,incThres=1)$N # enregistre le nbre d'individus
  Place <- records$trackRecord # enregistre les coordonnees x/y
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~ enregistrement distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dist <- as.matrix(list_data$particles$`Total movement`) # isole les donnees
  dist <- dist*4*10/125.6 # mise à l'echelle 100% et convertion pixels en mm (10mm = 125.6 pix)

  #~~~~~~~~~~~~~~~~~~~~~~~~~ donnes avec trajr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #initialisation des tableau pour enregistrer les donnees
  Tab_vit = matrix(ncol = 1, nrow = nbre) # vitesse
  Tab_sd = matrix(ncol = 1, nrow = nbre) # sd speed
  Tab_acc = matrix(ncol = 1, nrow = nbre) # acceleration
  Tab_arret = matrix(ncol = 1, nrow = nbre) # temps d'arret
  Tab_strait = matrix(ncol = 1, nrow = nbre) # stabilite (D/L)
  Tab_r = matrix(ncol = 1, nrow = nbre) # stabilite relative (r)
  Tab_sin = matrix(ncol = 1, nrow = nbre) # sinuosite
  Tab_Emax = matrix(ncol = 1, nrow = nbre) # Emax
  Tab_DC = matrix(ncol = 1, nrow = nbre) # DC
  Tab_SDDC = matrix(ncol = 1, nrow = nbre) # SDDC
  
  
  for (u in 1 : nbre){
    
    # Definition de x, y, et temp
    coords <- data.frame(x = Place[c(u),,c(1)], # x
                         y = Place[c(u),,c(2)], # y
                         time = seq((180/compteur), 180, by = 180/compteur)) # temps
    
    # création de la trajectoire 1 (pour : vitesse, sd vitesse, acceleration, temps d'arret)
    trj1 <- TrajFromCoords(coords, spatialUnits = "pixels", timeCol = 3)
    # Smooth before calculating derivatives
    smoothed <- TrajSmoothSG(trj1)
    
    # création de la trajectoire 2 (pour :stabilité (D/L), stabilite relative (r), sinuosite, Emax, DC, SDDC)
    trj2 <- lapply(1:1, function(i) TrajFromCoords(coords, spatialUnits = "pixels", timeCol = 3))
    # rediscretisation de la trajectoire
    reds <- lapply(1, function(ss) lapply(1:1, function(i) TrajRediscretize(trj2[[i]], ss)))
    
    #~~~~~~~~~~~~~~~~~~~~~ calcul vitesse et acceleration ~~~~~~~~~~~~~~~~~~~~~~

    # Calculate speed and acceleration
    derivs <- TrajDerivatives(smoothed)
    # sauvegarde les données
    Tab_vit[u] <- mean(derivs$speed*4*10/125.6) # vitesse
    Tab_sd[u] <- sd(derivs$speed*4*10/125.6) # vitesse dc
    Tab_acc[u] <- mean(derivs$acceleration*4*10/125.6) # acceleration
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~ calcul du temps d'arret ~~~~~~~~~~~~~~~~~~~~~~~~
    # Calcul des intervals
    intervals <- TrajSpeedIntervals(smoothed, slowerThan = 0.157) # en pix/s
    
    Tab_arret[u] <- sum(intervals$duration) # temps d'arret
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~ stabilite (D/L) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Tab_strait[u] <- sapply(reds, function(rtrjs) sapply(1:1, function(i) TrajStraightness(trj2[[i]]))) # D/L
    Tab_r[u] <- sapply(reds, function(rtrjs) sapply(1:1, function(i) Mod(TrajMeanVectorOfTurningAngles(trj2[[i]])))) # r

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sinuosite ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Tab_sin[u] <- sapply(reds, function(rtrjs) sapply(1:1, function(i) TrajSinuosity2(trj2[[i]]))) # sinuosite
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Emax ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Tab_Emax[u] <- sapply(reds, function(rtrjs) sapply(1:1, function(i) TrajEmax(trj2[[i]]))) # Emax
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Tab_DC[u] <- sapply(1:1, function(i) mean(TrajDirectionalChange(trj2[[i]])))
    Tab_SDDC[u] <- sapply(1:1, function(i) sd(TrajDirectionalChange(trj2[[i]])))
    
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~Compilation de toutes les donnees~~~~~~~~~~~~~~~~~~~~
  
  data <- cbind(conc[o],
                dist,
                Tab_vit,
                Tab_sd,
                Tab_acc,
                Tab_arret,
                Tab_strait,
                Tab_r,
                Tab_sin,
                Tab_Emax,
                Tab_DC,
                Tab_SDDC) # transforme les donnees en DF
  
  colnames(data) <- c("Concentration",
                      "distance [mm]",
                      "Vitesse moyenne [mm/s]",
                      "SD vitesse [mm/s]",
                      "accélération [mm/s^2]",
                      "Temps d'arret [s]", 
                      "straightness (D/L) [-]",
                      "alternative straigthness (r) [-]",
                      "sinuosite [-]",
                      "Emax [-]",
                      "DC [-]",
                      "SDDC [-]")
  
  data_all <- rbind(data_all,data) # rajoute les donnees aux tableaux
  
  if ( o == length(conc)){ 
    
    data_coords <- data_coords[,-1] # supprime la premiere colonne (NA)
    data_all <- data_all[-1,] # supprime la premiere ligne (NA)
    row.names(data_all) <- c(1:length(data_all[,1])) # change les index des lignes
    
    setwd("path")
    write.csv(data_all,file="donnees_comportements.csv")
    write.csv(data_coords,file="coordonnes.csv")
    setwd("path")
    
  }

}

#-------------------------------- fin de script --------------------------------
