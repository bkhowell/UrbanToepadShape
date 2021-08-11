#clear all current variables or functions from R memory and start (or restart) script/anayses
rm(list = ls(all = TRUE)) 

###### Howell Puerto Rican Anole Toe Pad Shape Preliminaruy Data Setup ######

## Necessry files for this R script

  # Howell PR anole shape.tps - located in the "KMW toe pad images" folder along side the images themselves
  # Howell PR anole semilandmark slider file.txt  - located in the "Final Analyses" folder

## Packages used in this script 

  #library(geomorph)
  packageVersion("geomorph")
  #3.3.1
  
  #library(shape) - nested within the included custom GM_IndvPlot() function 
  packageVersion("shape")
  #1.4.4

## Custom functions included in this script 

  # GM_curve_1st_last_remove() - removes the first and last semilandmarks from each individual's curves
  # GM_IndvPlot() - Returns plots of each specimen's landmarks and curves, color coated to facilitate the identificaiton of misplaced or misorderd landmarks or curves 

## My Current R and RStudio Versions 

  #check version of R
  getRversion()
  #R version 3.5.3

  #check R studio version
  require(rstudioapi)
  versionInfo()$version
  #R studio version 1.1.463

######
###### Load files ####

#Load geomoprh libary
library(geomorph)

#Set working directory
setwd("UrbanToepadShape")

#load TPS file
dat.raw<-readland.tps("Data/Howell PR anole shape.tps",specID="imageID", readcurve=T)
#the specID argument references what specimen names to use
#We didn't use the TPS file's "ID" variable, its values are arbitrary, please ignore
#The readcurve=T argument will treat landmarks in the TPS file under the CURVE headers as semiladnmarks

#TPS file should load in with no errors, noting "130 curve points detected per specimen and are appended to fixed landmarks." 
#This note represents our 13 curves per specimen with an initial 10 points per curve

#load slider file
sliders<-read.table("Data/Howell PR anole semilandmark slider file.txt")
#This 3 column tab delineated file dictates the relative relationship of each of our semilandmarks
#With each semilandmark in the center column, the left and rigth columns are the two other landmarks or semilandmarks the center semilandmark is connected to on either side




###### Remove First and Last Semilandmarks ####

#In order for the slider file to properly correspond with our semilandmarks, we must first remove the first and last semilandmark from each of our curves for all our individuals
#Each curve's first and last landmarks are redundant and overlap with the curve's anchoring landmarks

### MUST LOAD FUNCTION ###
#Run the GM_curve_1st_last_remove() function 

###### FUNCTION GM_curve_1st_last_remove() - removes 1st and last semilandmark points ####

## Arguments
# TPS_location - the file directory location and file name of the TPS file to be used, taking into account the working diretory 

GM_curve_1st_last_remove<-function(TPS_location){
  
  #To remove the first and last semilandmarks from each indivudals curves we need to read the TPS file as text, foregoing the use of geomorph's import function
  dat<-readLines(TPS_location)
  
  #We need to know number of landmarks, the number of curves, and  number of points in each curve
  #We can extract this data from the dat object
  #We assume all individuals have the same number of landmarks and same number of semilandmarks in their curves
  
  num_LM <- grep("LM=",dat) #Which lines have the text "LM" in them
  num_LM<-dat[num_LM[1]] #Get the first line with the text "LM" in it
  num_LM<-as.numeric(strsplit(num_LM, "=")[[1]][2]) #Split that line at the "=" and keep the text after the "="
  
  num_curves <- grep("CURVES=",dat) #Which lines have the text "CURVES" in them
  num_curves<-dat[num_curves[1]] #Get the first line with the text "CURVES" in it
  num_curves<-as.numeric(strsplit(num_curves, "=")[[1]][2]) #Split that line at the "=" and keep the text after the "="
  
  #Now that we know the number of curves, we need the number of points per curve 
  num_points<-matrix(ncol=1, nrow=num_curves) #make a blank matrix with a row for each curve
  
  for(c in 1:num_curves){ #Run a loop for each curve
    #c<-1
    out<-dat[grep("POINTS=",dat)[c]] #get the "POINTS = " line for the c-th curve
    num_points[c,1]<-as.numeric(strsplit(out, "=")[[1]][2]) #Split that line at the "=" and keep the text after the "="
  }
  
  #Now we need to figure out which landmarks to remove
  #for each curve, figure out the landmark numbers, which is the number of landmarks plus previous curves
  
  to_remove<-matrix(ncol=0, nrow=0) #make a blank matrix which will be populatioed by the semilandmarks we need to drop
  
  for(c in 1:num_curves){ #Run a loop for each curve
    #c<-1
    pts<-1:num_points[c,] #get a set of numbers for the c-th curve 1 through the number of semilandmarks in that curve
    pts<-pts+num_LM #add the number of landmarks
    pts<-pts+sum(num_points[0:(c-1),]) #add the number of semilandmarks in the previous curves
    out<-c(pts[1], pts[length(pts)]) #get the first and last semilandmark numbers of our focal curve
    to_remove<-c(to_remove, out) #add those 2 semilandmark numbers to our "to_remove" matrix
  }
  
  to_remove
  
}



to_remove<-GM_curve_1st_last_remove("Data/Howell PR anole shape.tps")

#We now have a list of the semilandmarks that we need to remove (to_remove) that is compatible with the way geomorph organises its data
dat.raw<-dat.raw[-to_remove,,] #drops first and last overlapping semilandmarks

#Check to make sure we now have the correct number of landmarks and semilandmarks and specimens
length(dat.raw[,1,1]) #number of total landmarks now
#Should be 123. 19 landmarks, and 8 semilandmarks per curve with 13 curves

length(dat.raw[1,1,]) #number of specimens
#should be 246 individuals



###### Geomorph Procrustes Alignment #####

#Procrustes alignment 
dat.super<-gpagen(dat.raw,curves=sliders, ProcD = F) #The ProcD=F argument denotes the use of bending energy to optimize the locaiton of our semilandmarks



###### Data Quality and Error Identification ####

plotOutliers(dat.super$coord) #Produce an outlier plot to check if any individuals need further attention

plotAllSpecimens(dat.super$coords) #Overlay all specimens using geomorph function

#You can produce a similar plot, overlaying all specimens, using native function that provide more plotting options 
plot(dat.super$coord[,1,], dat.super$coord[,2,], asp=1, col="grey", pch=19, cex=0.5, xlab=NA, ylab=NA, xaxt="n", yaxt="n")
points(dat.super$consensus, col="red", pch=19, cex=0.5)


#We use a custom fucntion to plot the unaligned and aligned landmarks and semilandmarks of each individual
### MUST LOAD FUNCTION ###
#Run the GM_IndvPlot() function

###### FUNCTION GM_IndvPlot () - plot/visualize landmarks and curves to assess order #####

## Arguments

# gpagen_output - This fucntion uses geomorphs gpagen function output (aligned data)
# TPS_raw - It also unaligned data imported via geomorph (after first and last semilandmarks have been removed)
# slider_dat - It also uses the same landmark slider file
# plot_destination - The user need to dictate a save location for the created plots to go 

GM_IndvPlot<-function(gpagen_output, TPS_raw, slider_dat, plot_destination){
  
  library(shape) #This function requires the "Arrows" function from the shape library
  #gpagen_output<-dat.super
  #TPS_raw<-dat.raw
  #slider_dat<- sliders
  #plot_destination<-"Landmark Diagnostics/"
  
  
  curves<-as.list(as.data.frame(t(slider_dat))) #use the slider file to get the semilandmark relationships for each curve 
  
  #we need to build a set of landmarks and semilandmarks for each curve
  for(r in 1:length(curves)){ #Run a loop for each cruve
    #r<-2
    second_tolast<-lapply(curves, function(x) x[length(x) -1]) #get each semilandmark from the slider file (ie the 2nd to last elements of each object in the list)
    
    if(curves[[r]][1] %in% second_tolast){       #for each semilandmark, if the point before it is also a semilandmark... 
      togo<-which(second_tolast==curves[[r]][1]) #Find the preceding semilandmark in the slider file
      curves[[togo]]<-c(curves[[togo]], curves[[r]][3]) #ammend the landmark (or semilandmark) that follows the current landmark to the preceding semilandmark's slider line
      curves[[r]]<-NA #clear the current semilandmarks line in the slider file
    }
  }
  
  curves<-curves[!is.na(curves)] #remove the emply lines from the slider file
  num_curves<-length(curves) #get the number of curves
  names(curves)<-paste("semi", 1:num_curves, sep="") #rename the objects sequentially
  #This should result in a list of objects, one for each curve, that lists the landmarks and semilandmarks in that curve in the correct order
  
  
  hard<-c(paste(lapply(curves, function(x) x[1])), paste(lapply(curves, function(x) x[length(x)]))) #get the first and last landmark from each curve
  hard<-sort(unique(as.numeric(hard))) #transform them into numeric variables and sort them
  
  #start translating variables 
  gpagen_output<-gpagen_output$coords #extract just the coordinate data from the gpagen output data
  nsamples<-dim(TPS_raw)[3] #get the number of individuals that we'll need plots for
  
  for (spec_num in 1:nsamples) {  #start cycling through individuals
    #spec_num<-1
    
    #Get that individual's data
    specimen_aligned<-gpagen_output[,,spec_num]
    specimen_raw<-TPS_raw[,,spec_num]
    
    Indv<-dimnames(gpagen_output)[[3]][spec_num] #get the individuals name 
    name1<-Indv
    #These sample names cannot contain pathway information. They must be cropped to only inlucde the specimen name
    #geomorph automatically removes the filename extension from the specimen names 
    
    
    #set raw all positive, near 0
    specimen_raw[,1]<-specimen_raw[,1]-(min(specimen_raw[,1])) #shift x values
    specimen_raw[,2]<-specimen_raw[,2]-(min(specimen_raw[,2])) #shift y values
    
    
    #start making plots
    pdf(file=(file.path(paste(plot_destination, name1, ".pdf", sep=""))), width=7, height = 4)
    
    
    par(mfrow=c(1,2), mar=c(1,3,1,1), pty="s") 
    
    
    #unaligned points first 
    plot(specimen_raw[,1], specimen_raw[,2], xlim=c(min(specimen_raw), max(specimen_raw)), ylim=c(min(specimen_raw), max(specimen_raw)), xlab=NA, ylab=NA, col="white", main="Unaligned",  cex.axis=0.5)
    mtext("(mm)", side=1, line=2) #this assumes your scale in TPSDig was in millimeters
    mtext("(mm)", side=2, line=2) #this assumes your scale in TPSDig was in millimeters
    
    
    mtext(text=Indv, cex=0.7, side = 3, outer = T, line = -2)
    
    for (l in 1:num_curves) { #plot each curve
      #l<-1
      lines(x=c(specimen_raw[curves[[l]],1]), y=specimen_raw[curves[[l]],2], col=rainbow(num_curves)[l], type="o", pch=19, cex=0.3)
      Arrows(x0=specimen_raw[curves[[l]][length(curves[[l]])-1],1], x1=specimen_raw[curves[[l]][length(curves[[l]])],1], y0=specimen_raw[curves[[l]][length(curves[[l]])-1],2], y1=specimen_raw[curves[[l]][length(curves[[l]])],2], col=rainbow(num_curves)[l], arr.type="triangle", arr.length = 0.15, arr.width=0.15, segment=F)
    }
    #This function plots each curve, using a unique sequencial color for each curve, ending with an arrow
    
    points(specimen_raw[hard,], pch=19, cex=0.5) #add landmarks
    text(specimen_raw[hard,], label=hard, cex=0.5, pos=3, offset=0.2) #adds  labels to the landmarks
    
    
    
    #aligned data
    plot(specimen_aligned[,1], specimen_aligned[,2], xlim=c(min(specimen_aligned), max(specimen_aligned)), ylim=c(min(specimen_aligned), max(specimen_aligned)), xlab=NA, col="white", main="Aligned",  cex.axis=0.5)
    
    for (l in 1:num_curves) { 
      #l<-1
      lines(x=c(specimen_aligned[curves[[l]],1]), y=specimen_aligned[curves[[l]],2], col=rainbow(num_curves)[l], type="o", pch=19, cex=0.3)
      Arrows(x0=specimen_aligned[curves[[l]][length(curves[[l]])-1],1], x1=specimen_aligned[curves[[l]][length(curves[[l]])],1], y0=specimen_aligned[curves[[l]][length(curves[[l]])-1],2], y1=specimen_aligned[curves[[l]][length(curves[[l]])],2], col=rainbow(num_curves)[l], arr.type="triangle", arr.length = 0.15, arr.width=0.15, segment=F)
    }
    #This function plots each curve, using a unique sequencial color for each curve, ending with an arrow
    
    points(specimen_aligned[hard,], pch=19, cex=0.5)
    text(specimen_aligned[hard,], label=hard, cex=0.5, pos=3, offset=0.2) #adds labels to the hard landmarks
    
    dev.off()
    
  }
  
}

#This function produces a multiplot .pdf file for each individual. You must have a destination folder already in place to receive these .pdf plots
#I use the folder "Plots" within our "Data" folder to receive these plots  

#Run the function to visualize each individuals landmarks and curves
GM_IndvPlot(gpagen_output = dat.super, TPS_raw = dat.raw, slider_dat = sliders, plot_destination = "Data/Plots/")

#This function allows users to quickly diagnose individuals in which the landmarks or curves may be out of order
#The function uniquely colors each curve highlighting the order in which the curves should appear. Curves also include an arrow highlighting the order of points
#Since each curve is anchored between two landmarks, if any curves or anchoring landmrks are out of order, this will be evident by misplaced lines connecting them

###### Import Field Data ####

##Our field data includes each specimen's:

#Unique specimen name 
#The municipality it was collected  
#The site in which that it was collected
#The habitat category of the site in which it was collected (urban or forest)
#Body length (SVL) in millimeters
#Toe pad area in square millimeters. 
#Number of toe pad lamellae 

#The toes used for to measure area and count lamellae were the same longest rear toes we used for our geometric morphometrics

#Import field data
dat.field<-read.delim("Data/Howell Field Data.txt", as.is = T, header=T) #as.is=T to prevent R from automatically assigning text variable as factors, header=T to automatically assign column names

#Confirm that the specimen names in our landmark data and field data match and are in the same order
dimnames(dat.super$coords)[[3]]==dat.field$Specimen
#should all be true


###### Export Data ####

#Combine our pre-alignment shape data (dat.raw), aligned shape data (dat.super), and our field data (field.dat) into one R object for export

Howell_dat<-list(dat.raw, dat.super, dat.field)
names(Howell_dat)<-c("dat.raw", "dat.super", "dat.field") #assign names to each inlcuded datasets 

save(Howell_dat, file="Howell R data.Rdata")







