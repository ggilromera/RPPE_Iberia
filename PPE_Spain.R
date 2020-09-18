library(DEoptim)
setwd("~/IPE/PragaWork/EMPDIFN/5.PPE")
dwm="lsm unstable"


fname_lc <- "vegetation.xlsx"                                                         # plant cover (lc stands for landcover)

fname_pollen <- "pollen.xlsx"                                               # pollen data

fname_parameter <- "flsp_param.xlsx"                                     # parameter

outfile <- "results.csv"                                                        # results file name  



#################################################

# load data

#################################################



library(readxl)

sheets <- readxl::excel_sheets(fname_lc)                                      # read land (plant) cover

lc <-    lapply(sheets, function(X) readxl::read_excel(fname_lc, sheet = X))

names(lc) <- sheets

pollen <- readxl::read_excel(fname_pollen, sheet = 1, col_names = TRUE)       # read pollen data

parameter <- readxl::read_excel(fname_parameter, sheet = 1, col_names = TRUE) # read further parameters

pollen[is.na(pollen)] <- 0
lc[is.na(lc)] <- 0


######################################

# prepare parameters 

######################################



vg <- parameter[,"fallspeed"]                                       # read fall speed

nTaxa <- nrow(pollen)                                               # read number of taxa



rad <- lc[[1]][,1]                                                  # read vector of radii

radii <- c(0.5, rad$Radius)                                                # add zero for first ring

areaCircles <- pi * (radii^2); areaRings <- abs(diff(areaCircles))  # calculate ring areas  



# hier kann Radius raus...

lcPer <- lapply(lc, function(x) x/areaRings)                        # proportions of landcover 

ri <- lapply(vg$fallspeed, FUN=ringInflux, dwm=dwm, radii=radii)              # species specific weighting factor for each ring







######################################

# core function

# for repeated calculations of PPEs

# not yet complete

######################################



n <- length(ri)

lcWeighted <- lapply(1:n,function(i,m,v){m[[i]]*v[[i]]},m=lcPer,v=ri)      # multiply landcover of each taxon with corresponding weigting factor



lower <- rep(0.1, nTaxa); upper <- rep(20, nTaxa)                   # set boundaries for PPEs in optimization



# add random error to pollen data (by default: random drawing with rnorm from the given composition)

pollenfun=rmultinom_reveals


## loop for repetitions with subset of vegetation data

for (s in 1:5){
  
  
  
  # add error to pollen data
  
  pollenPer <- apply(pollen[-1], 2, pollenfun, n=1)
  
  write.table(s, file="pcMM_Results_increasing radii_check.csv",col.names=F, append = TRUE, sep=";", dec=",")
  
  
  
  
  
  # loop over radii (the number was 53 in Tasmanian data) 35
  
  for (i in 2:8){
    
    ma <- list()
    
    for (j in 1:nTaxa){
      
      ma[[j]] <- as.matrix(lcWeighted[[j]][1:i,])
      
    }
    
    
    
    # calls DeOptim
    
    oX <- DEoptim(optimPPE.c, lower, upper, ma=ma, pollenPer=pollenPer, nTaxa=nTaxa,
                  
                  DEoptim.control(strategy=2, itermax = 20000, reltol=0.01, steptol=500, VTR = 0, parallelType=1))
    
    
    
    write.table(round(t(c(oX$optim$bestval,oX$optim$bestmem)), digits=3), file=outfile,col.names=F, append = TRUE, sep=";", dec=",")
    
  }
  
  
  
}






save(oX, file="oX.RData")
# now follow the three functions needed, the later 2 are already part of the Disqover package





##################################################################################

#####    calculates distance between modeled and empiric pollen percentages

##################################################################################





#optimPPE = function(ppe, lcWeighted, pollenPer, nTaxa){

optimPPE = function(ppe, ma, pollenPer, nTaxa){
  
  

  # dies hier bremst - obwohl gar nicht viel passiert
  
  lcWeightedPPE <- lapply(1:nTaxa,function(i,m,v){m[[i]]*v[i]},m=ma,v=ppe)  
  
  
  
  pollenInflux <- lapply(lcWeightedPPE, function(x) colSums(x)) 
  
  
  
  influx <- do.call(rbind, pollenInflux) 
  
  influx <- influx[,-1]
  
  
  
  pollenMod <- 100*prop.table(influx, 2)
  
  
  
  #calculate distance between modelled and empiric pollen data
  
  #dist <- sum(((pollenMod-pollenPer)^2)/((pollenPer+1)*ppe))    # scheint hier nicht gut zu laufen, keine Idee warum
  
  dist <- sum(((pollenMod-pollenPer)^2)/((pollenPer+1))) 
  
  return(dist)
  
}





# to speed things up --> compile function in c

library(compiler)

optimPPE.c <- cmpfun(optimPPE)





##################################################################################

##### distance weighting factor for each ring for one taxon

##################################################################################



# this Version is differenf from the package!! and possibly better





ringInflux = function(vg, dwm, radii){
  
  
  
  ###
  
  # die Abfrage muss nur einmal passieren, nicht jedes mal!
  
  ###
  
  
  
  # branch with dispersal model (for parameters of the GPM see Jackson & Lyford 1999)
  
  if (dwm=="gpm neutral") {
    
    #disModel<-GPM(fallSpeed=vg, cz=0.12, n=0.25, u=3)
    
    a <- lapply(radii, GPM, fallSpeed=vg, cz=0.12, n=0.25, u=6.5)
    
    airborne <- do.call(rbind, a)
    
    
    
  } else if (dwm=="gpm unstable") {
    
    #disModel<-GPM(fallSpeed=vg, cz=0.21, n=0.20, u=3)
    
    a <- lapply(radii, GPM, fallSpeed=vg, cz=0.21, n=0.20, u=3)
    
    airborne <- do.call(rbind, a)
    
    
    
  } else if (dwm=="1overd") {
    
    # nicht ausgearbeitet, oder?
    
    disModel<-OneOverD.f()
    
    airborne <- predict(disModel, radii)
    
    
    
  } else if (dwm=="lsm unstable")   {
    
    load("Kuparinen unstable airborne2.rda") # load model, which is prepared using a loess smoother
    
    
    
    vgr<-100*round(vg,2)+4 #find right column
    
    if (vg==0.002) vgr <- 3
    
    if (vg==0.001) vgr <- 2
    
    
    
    disModel<-kup_unstable_airborne[[vgr]]
    
    airborne <- predict(disModel, radii)   # pollen airborne at each circle 
    
    
    
    #} else if (dwm=="alternative")   { # with the following option, alternative models may be used...
    
    #  load("alternative.rda") # load model, which is prepared using a loess smoother
    
    #  vgr<-100*round(vg,2)+1 #find right column
    
    #  disModel<-alternative[[vgr]]
    
  } else stop("no valid distance weighting method selected")
  
  
  
  ringInflux <- abs(diff(airborne)) 
  
  
  
  return(ringInflux)
  
} 





##########################################################################

##### Gaussian plume model (GPM) - Sutton model follwing Prentice (1985) 

##########################################################################



#  a <- GPM(fallSpeed=vg, cz=0.12, n=0.25, u=3)





GPM = function (x, fallSpeed, cz, n, u){
  
  y <- exp((-4*fallSpeed*x^(n/2))/(n*u*sqrt(pi)*cz))
  
  return(y)
  
}







##################################################################################

#####    error in pollen data

##################################################################################





## this function is already part of the package...., just needed here for practical reasons



rmultinom_reveals<-function(n=1,pollen,...)
  
{
  
  100*rmultinom(n, sum(pollen), pollen/sum(pollen))/sum(pollen)
  
}




