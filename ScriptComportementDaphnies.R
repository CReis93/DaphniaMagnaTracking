###############################
####
####      Analyse du comportement_Exp daphnies
####      Auteur-trice : Christophe Reis & Floriane Tisserand
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
setwd("/work/FAC/FGSE/IDYST/nchevre/comp_daph/")

#~~~~~~~~~~~~~~~~~~~~ selection de l'experience ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ne runner que la ligne correspondantea l'experience
#exp = "low_conc"
exp = "high_conc_24h"
#exp = "high_conc_48h"

#~~~~~~~~~~~~~~~~~~ informations vidéos et DL vidéos ~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(exp == "low_conc"){ # videos low conc.

  conc <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10") # concentrations
  individus <- c(14,14,8,7,6,9,8,9,14,8) # nbre d'individus
  thres <- c(-0.115,-0.115,-0.115,-0.115,-0.115,-0.115,-0.115,-0.115,-0.115,-0.18) # seuil
  images.manq <- c(0,0,0,0,0,0,0,0,0,0) # vct pour image manquante
  # telecharge les videos low conc.
  list_mov <- list.files(path = "videos/brute", pattern = (".mp4"), full.names = TRUE)
  
}

if(exp == "high_conc_24h"){ # videos high conc. 24h
  
  conc <- c("c0","c1","c2","c3","c4") # concentrations
  individus <- c(8,9,9,10,7) # nbre d'individus
  thres <- c(-0.115,-0.115,-0.115,-0.115,-0.115) # seuil
  images.manq <- c(0,0,0,0,0) # vct pour image manquante
  # telecharge les videos high conc. 24h
  list_mov <- list.files(path = "videos/brute/24H", pattern = (".mp4"), full.names = TRUE)
  
}

if(exp == "high_conc_48h"){ # videos high conc. 48h
  
  conc <- c("c0","c1","c2","c3") # concentrations
  individus <- c(5,10,10,10) # nbre d'individus
  thres <- c(-0.115,-0.115,-0.115,-0.115) # seuil
  images.manq <- c(0,0,0,0) # vct pour image manquante
  # telecharge les videos high conc. 48h 48H
  list_mov <- list.files(path = "videos/brute/48H", pattern = (".mp4"), full.names = TRUE)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ variables generales ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# nombre de FPS
FramePS <- 15

# creation tableau all_data pour sauvegarder les donnees
data_all = matrix(ncol = 4)
colnames(data_all) <- c("Concentration",
                        "Vitesse moyenne [mm/s]",
                        "Temps d'arret [s]", 
                        "sinuosite")

#~~~~~~~~~~~~~~~~~~~~~~~~~ remise à 0 du tableau des données ~~~~~~~~~~~~~~~~~~

file.remove(list.files(path = "resultats/imagesChemin", 
                       pattern = (".pdf"), 
                       full.names = TRUE))
file.remove(list.files(path = "resultats/videoChemin", 
                       pattern = (".png"), 
                       full.names = TRUE))
file.remove(list.files(path = "resultats/", 
                       pattern = (".csv"), 
                       full.names = TRUE))

#~~~~~~~~~~~~~~~~~~ Convert video to images sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (o in 1:length(conc)){ # pour les videos de Floriane

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
  
  resume <- summary(partIden) # sauvegarde le resultat de identifyParticules

  miss <- which(resume$n != individus[o]) # prend les images avec ind. manquant
  images.manq[o] <- length(miss) # sauvegarde le nbre d'images
  N.images <- N.images[-c(miss)] # enleves les images avec des individus en trop
  
  if(resume$particles[2] != 0){

    # On recommence tout sans les images "defectueuses" ........................
  
    allFullImages <- loadImages(dirPictures = "videos/images",
                              nImages = N.images)
    stillBack <- createBackground(allFullImages,method="mean")
    allImages <- subtractBackground(stillBack)
    partIden <- identifyParticles(allImages,
                                threshold = thres[o],
                                pixelRange = c(5,500)) # reporter la valeur seuile
    print(summary(partIden)) # montre le resultat final pour controle
  
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~ Enregistrement des mouvements ~~~~~~~~~~~~~~~~~~~~~

  records <- trackParticles(partIden,L=200,R=3) # traque et lie les points. 

  #~~~~~~~~~~~~~~~~~~~ créer l'image de suivit et la video ~~~~~~~~~~~~~~~~~~~~~
  
  setwd("/work/FAC/FGSE/IDYST/nchevre/comp_daph/resultats/imagesChemin")
  ext <- ".pdf"
  name <- "imageSuivis_C"
  num <- as.character(o)
  fichier <- paste(name,num,ext,sep="")
  pdf(file = fichier)
  par(mfrow=c(1,1))
  plot(records,type="trajectories",incThres=1) # affichage du resultat
  dev.off()
  setwd("/work/FAC/FGSE/IDYST/nchevre/comp_daph/")

  #plot(records,type="animation",incThres=1,
       #path="resultats/videoChemin/",
       #libavpath = "libav/usr/bin/avconv.exe") # creation video de suivis
  
  #name <- as.character(o) # trasnforme le numero de la video en caractere
  #ext <- ".mp4" # creation de l'extension
  #chemin <- "resultats/videoChemin/" # creation du chemin
  #fichier <- paste(chemin,name,ext,sep="") # merge des trois string
  
  #file.rename(from="resultats/videoChemin/animation.mp4", to = fichier) # renome le fichier
  
  #~~~~~~~~~~~~~~~~~~~~~~~ Enregistrement des donnees ~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  resultat <- summary(records, incThres = 1)
  list_data <- Map(as.data.frame, resultat) # extrait les donnees
  Move <- as.data.frame(list_data$particles$`Total movement`) # isole les donnees
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ converstion en vitesse ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Move <- Move*4 # mise à l'echelle 100% 
  Move <- Move*10/125.6 # conversion pixels en mm (10mm = 125.6 pix)
  Vit <- Move/180 # conversion en vitesse (mm/s) avec video de 3min (180 sec.)
  
  #~~~~~~ donnees utile pour le temps d'arret et sinuosite ~~~~~~~~~~~~~~~~~~~~~
  
  nbre <- summary(records,incThres=1)$N # enregistre le nbre d'individus
  step <- compteur-length(miss) # enregistre le nombre d'images (steps)
  Place <- records$trackRecord # enregistre les coordonnees x/y
  
  # creation des tableau pour sauvegarde des resultats
  Tab_arret = matrix(ncol = 1, nrow = nbre) # creation d'une matrice temps d'arret
  Tab_sin = matrix(ncol = 1, nrow = nbre) # creation d'une matrice sinuosite
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ calcul du temps d'arret ~~~~~~~~~~~~~~~~~~~~~~~~~~

  for (l in 1: nbre){
    
    # Definition de x, y, et temp
    
    coords <- data.frame(x = Place[c(l),,c(1)], # x
                         y = Place[c(l),,c(2)], # y
                         time = seq((180/step), 180, by = 180/step))# temps
    
    # mise a l'echelle de pixel en metre
    trj <- TrajFromCoords(coords, spatialUnits = "pixels")

    #  10 mm (=   1 cm) =  125.6 pixels dans la video
    trj <- TrajScale(trj, 10 / 125.6, "mm")

    # analyse de la trajectoire
    
    # Calcul des intervals
    intervals <- TrajSpeedIntervals(trj, slowerThan = 0.05) # en mm/s
    
    #plot(intervals)
    Tab_arret[l] <- sum(intervals$duration)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sinuosite ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  for (m in 1: nbre){
    
    # Definition de x, y, et temp
    trajectoire <- data.frame(x = Place[c(m),,c(1)], # x
                            y = Place[c(m),,c(2)], # y
                            time = seq((180/step), 180, by = 180/step)) # temps
  
    # mise a l'echelle de pixel en metre
    trajectoire2 <- TrajFromCoords(trajectoire, spatialUnits = "pixels")
  
    # Calcule le changement de direction
    Tab_sin[m] <- TrajSinuosity2(trajectoire2)
  
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~Compilation de toutes les donnees~~~~~~~~~~~~~~~~~~~~
  
  data <- cbind(conc[o],Vit,Tab_arret,Tab_sin) # transforme les donnees en DF
  colnames(data) <- c("Concentration","Vitesse moyenne [mm/s]","Temps d'arret [s]", "sinuosite")
  data_all <- rbind(data_all,data) # rajoute les donnees aux tableaux
  
  if ( o == length(conc)){ 
    
    data_all <- data_all[-1,] # supprime la premiere ligne (NA)
    row.names(data_all) <- c(1:length(data_all[,1])) # change les index des lignes
    
    # Sortie images manquante et pourcentage
    images.manq.pour <- images.manq/compteur*100
    donne.manq <- cbind(images.manq,images.manq.pour)
    
    setwd("/work/FAC/FGSE/IDYST/nchevre/comp_daph/resultats")
    write.csv(data_all,file="donnees_comportements.csv")
    write.csv(donne.manq,file="imagesManquante.csv")
    write.csv(seuil_vitesse,file="seuil_vitesse.csv")
    setwd("/work/FAC/FGSE/IDYST/nchevre/comp_daph/")
    
  }

}

#-------------------------------- fin de script --------------------------------
