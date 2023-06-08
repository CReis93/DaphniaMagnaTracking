###############################
####
####      Behaviour analysis - Daphnia Magna expsed to pesticide
####      Authors : Christophe Reis & Floriane Tisserand
####      Version : 19/04/2023
####
####
####      sources : 
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

# define the WD
setwd("path") # path == working directory on computer

#~~~~~~~~~~~~~~~~~~~~~~~ Experience selection ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# run only one line out of the three
exp = "low_conc"
#exp = "high_conc_24h"
#exp = "high_conc_48h"

#~~~~~~~~~~~~~ informations about the videos et DL videos ~~~~~~~~~~~~~~~~~~~~~~

if(exp == "low_conc"){ # videos low conc.
  
  conc <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10") # concentrations
  individus <- c(12,14,7,7,6,9,8,9,13,8) # number of individuals
  thres <- c(-0.12,-0.18,-0.16,-0.18,-0.21,-0.21,-0.115,-0.14,-0.18,-0.27) # threshold
  list_mov <- list.files(path = "videos/brute/", pattern = (".mp4"), full.names = TRUE) # download (DL) the videos
  
}

if(exp == "high_conc_24h"){ # videos high conc. 24h
  
  conc <- c("c0","c1","c2","c3","c4") # concentrations
  individus <- c(8,9,9,10,7) # number of individuals
  thres <- c(-0.115,-0.08,-0.115,-0.3,-0.3) # threshold
  list_mov <- list.files(path = "videos/brute/24H", pattern = (".mp4"), full.names = TRUE) # DL the videos
  
}

if(exp == "high_conc_48h"){ # videos high conc. 48h
  
  conc <- c("c0","c1","c2","c3") # concentrations
  individus <- c(5,9,10,10) # number of individuals
  thres <- c(-0.115,-0.115,-0.07,-0.115) # threshold
  list_mov <- list.files(path = "videos/brute/48H", pattern = (".mp4"), full.names = TRUE) # DL the videos
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ general variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FPS number
FramePS <- 15

# creation of results storage matrices
data_all = matrix(ncol = 12)
colnames(data_all) <- c("Concentration", "distance [mm]", "Vitesse moyenne [mm/s]", "SD vitesse [mm/s]",
                        "accélération [mm/s^2]", "Temps d'arret [s]", "straightness (D/L) [-]",
                        "alternative straigthness (r) [-]", "sinuosite [-]", "Emax [-]", "DC [-]",
                        "SDDC [-]")

#~~~~~~~~~~~~~~~~~~~~~~~~ clean all the directory with result ~~~~~~~~~~~~~~~~~

file.remove(list.files(path = "resultats/imagesChemin", pattern = (".pdf"), full.names = TRUE))
file.remove(list.files(path = "resultats/videoChemin", pattern = (".mp4"), full.names = TRUE))
file.remove(list.files(path = "resultats/", pattern = (".csv"), full.names = TRUE))

#~~~~~~~~~~~~~~~~~~ Convert video to images sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (o in 1:length(conc)){
 
  print(paste("On commence la concentration num :",o)) # checking the start of a new loop
  
  # deletes old images from the directory
  file.remove(list.files(path = "videos/images", 
                         pattern = (".png"), 
                         full.names = TRUE))
  
  # cut to the video into images
  av_video_images(list_mov[o], 
                  destdir = "videos/images", 
                  format = "png", fps = FramePS)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DL the images ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  compteur <- length(list.files("videos/images"))-1 # cunt the number of images
  N.images <- c(1:compteur) # define a vector containing the number of images
  allFullImages <- loadImages (dirPictures = "videos/images",
                               nImages = N.images)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ background noise suppression ~~~~~~~~~~~~~~~~~~~~~
  
  #1) creation of background noise image
  stillBack <- createBackground(allFullImages,method="mean")
  
  #2) erase the background noise image
  allImages <- subtractBackground(stillBack)
  
  #~~~~~~~~~~~~~~~~~~~ identification of mobile particles ~~~~~~~~~~~~~~~~~~~~~
  
  partIden <- identifyParticles(allImages,
                                threshold = thres[o],
                                pixelRange = c(5,500))
  
  print(summary(partIden)) # show the final result for controle
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~ recording movements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  records <- trackParticles(partIden,L=200,R=3) # tracks down and connects the dots. 
  
  #~~~~~~~~~~~~~~~~~~~~ create tracking image and video ~~~~~~~~~~~~~~~~~~~~~~~
  
  setwd("path/resultats/imagesChemin")
  pdf(file = paste("imageSuivis_C",as.character(o),".pdf",sep=""))
  par(mfrow=c(1,1))
  plot(records,type="trajectories",incThres=1) # plot the result
  dev.off()
  setwd("path")
  
  plot(records,type="animation",incThres=1,
  path="resultats/videoChemin/",
  libavpath = "libav/usr/bin/avconv.exe") # creation of the animation
  file.rename(from="resultats/videoChemin/animation.mp4", 
              to = paste("resultats/videoChemin/",as.character(o),".mp4",sep="")) # rename the file
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~ Registration of coordinates ~~~~~~~~~~~~~~~~~~~~~~
  
  coordonnes <- matrix(ncol = (individus[o]*2), nrow = (length(records$trackRecord[1,,1]))) # coordinate matrix creation
 
  wu <- 1 # temporary variables
  wuu <- 1
  
  for(i in 1 : (individus[o]*2)){ # save the coordinates
    
    if(wuu==3){wu <- wu + 1}
    if(wuu==3){wuu <- 1}
    
    coordonnes[,i] <- records$trackRecord[wu,,wuu]
    
    wuu <- wuu + 1
    
  }
  
  if(o==1){data_coords = matrix(nrow = compteur)} # create ending table for coordinates
  data_coords <- cbind(data_coords,coordonnes) # add the new coordinates
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ save the data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  resultat <- summary(records, incThres = 1)
  list_data <- Map(as.data.frame, resultat) # extract the data
  
  #~~~~~~~~~~~~~~~~~~ data useful for behaviour analysis ~~~~~~~~~~~~~~~~~~~~~~~
  
  nbre <- summary(records,incThres=1)$N # save the number of individuals
  Place <- records$trackRecord # save the coordinates X and Y
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~ distance recording ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  dist <- as.matrix(list_data$particles$`Total movement`) # isolates data
  dist <- dist*4*10/125.6 # scale to 100% and convert pixels to mm (10mm = 125.6 pix)

  #~~~~~~~~~~~~~~~~~~~~~~~~~ data acquisition with trajr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # initializing tables to save data
  Tab_vit = matrix(ncol = 1, nrow = nbre) # speed
  Tab_sd = matrix(ncol = 1, nrow = nbre) # standard-deviation of speed
  Tab_acc = matrix(ncol = 1, nrow = nbre) # acceleration
  Tab_arret = matrix(ncol = 1, nrow = nbre) # stopping time
  Tab_strait = matrix(ncol = 1, nrow = nbre) # stability (D/L)
  Tab_r = matrix(ncol = 1, nrow = nbre) # relative stability (r)
  Tab_sin = matrix(ncol = 1, nrow = nbre) # sinuosity
  Tab_Emax = matrix(ncol = 1, nrow = nbre) # Emax
  Tab_DC = matrix(ncol = 1, nrow = nbre) # DC
  Tab_SDDC = matrix(ncol = 1, nrow = nbre) # SDDC
  
  for (u in 1 : nbre){
    
    # definition of x, y, and time
    coords <- data.frame(x = Place[c(u),,c(1)], # x
                         y = Place[c(u),,c(2)], # y
                         time = seq((180/compteur), 180, by = 180/compteur)) # time
    
    # create trajectory 1 (for : speed, sd-speed, acceleration, stop time)
    trj1 <- TrajFromCoords(coords, spatialUnits = "pixels", timeCol = 3)
    # Smooth before calculating derivatives
    smoothed <- TrajSmoothSG(trj1)
    
    # create trajectory 2 (for :stability (D/L), relative stability (r), sinuosity, Emax, DC, SDDC)
    trj2 <- lapply(1:1, function(i) TrajFromCoords(coords, spatialUnits = "pixels", timeCol = 3))
    # trajectory's rediscretisation
    reds <- lapply(1, function(ss) lapply(1:1, function(i) TrajRediscretize(trj2[[i]], ss)))
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~ speed ans acceleration ~~~~~~~~~~~~~~~~~~~~~~~~~

    # calculate speed and acceleration
    derivs <- TrajDerivatives(smoothed)
    # save data
    Tab_vit[u] <- mean(derivs$speed*4*10/125.6) # speed
    Tab_sd[u] <- sd(derivs$speed*4*10/125.6) # sd-speed
    Tab_acc[u] <- mean(derivs$acceleration*4*10/125.6) # acceleration
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ stopping time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # calculating intervals
    intervals <- TrajSpeedIntervals(smoothed, slowerThan = 0.157) # speed in [pix/s]
    
    Tab_arret[u] <- sum(intervals$duration) # stopping time
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~ stability (D/L) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Tab_strait[u] <- sapply(reds, function(rtrjs) sapply(1:1, function(i) TrajStraightness(trj2[[i]]))) # D/L
    Tab_r[u] <- sapply(reds, function(rtrjs) sapply(1:1, function(i) Mod(TrajMeanVectorOfTurningAngles(trj2[[i]])))) # r

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sinuosity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Tab_sin[u] <- sapply(reds, function(rtrjs) sapply(1:1, function(i) TrajSinuosity2(trj2[[i]]))) # sinuosity
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Emax ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Tab_Emax[u] <- sapply(reds, function(rtrjs) sapply(1:1, function(i) TrajEmax(trj2[[i]]))) # Emax
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Tab_DC[u] <- sapply(1:1, function(i) mean(TrajDirectionalChange(trj2[[i]])))
    Tab_SDDC[u] <- sapply(1:1, function(i) sd(TrajDirectionalChange(trj2[[i]])))
    
  }
    
  #~~~~~~~~~~~~~~~~~~~~~~~~ Compilation of all data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  data <- cbind(conc[o], dist, Tab_vit, Tab_sd, Tab_acc, Tab_arret, Tab_strait,
                Tab_r, Tab_sin, Tab_Emax, Tab_DC, Tab_SDDC) # transforme data in dataframe
  
  colnames(data) <- c("Concentration", "distance [mm]", "Vitesse moyenne [mm/s]", "SD vitesse [mm/s]",
                      "accélération [mm/s^2]", "Temps d'arret [s]", "straightness (D/L) [-]",
                      "alternative straigthness (r) [-]", "sinuosite [-]", "Emax [-]", "DC [-]",
                      "SDDC [-]")
  
  data_all <- rbind(data_all,data) # add data to the final table
  
  if ( o == length(conc)){ # run only when all the videos are analyzed
    
    data_coords <- data_coords[,-1] # erase the first line (NA)
    data_all <- data_all[-1,] # erase the first line (NA)
    row.names(data_all) <- c(1:length(data_all[,1])) # change the index of the table
    
    setwd("path/resultats") # save all results in a .csv file
    write.csv(data_all,file="donnees_comportements.csv")
    write.csv(data_coords,file="coordonnes.csv")
    setwd("path")

  }
  
}

#-------------------------------- end of script --------------------------------

