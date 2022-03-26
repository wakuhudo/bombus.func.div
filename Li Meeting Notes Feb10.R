library(tidyverse)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(vegan)
library(ggvegan)
library(gstat)
library(lattice)
library(nlme)
library(shape)
library(bipartite)
library(labdsv)
library(Matrix)
library(rgl)
library(ggpubr)
library(stringr)
library(mctest)
library(geomorph)
library(alphahull)
library(sf)

#read wing data
temp <- list.files("/Users/Apple/Desktop/left_2_13_15LM" ,pattern="*.TPS", full.names = TRUE)
tempr<- list.files("/Users/Apple/Desktop/right_2_13_15LM" ,pattern="*.TPS", full.names = TRUE)
wingdata<-readmulti.tps(temp, specID=c("imageID"))
wingdata
wingdatar[,1,] <- wingdatar[,1,] * -1 
bees_R_2D <- two.d.array(wingdatar)
bee_combined <- rbind(bees_L_2D, wingdatar) %>% arrayspecs(p=13, k=2) 
bee_gpa1 <- gpagen(bee_combined)
plot(bee_gpa1)
#bee_gpa <- gpagen(wingdata,ProcD = T)
#csize <- bee_gpa$Csize %>% as.matrix()
#csize <- data.frame(ImageID = row.names(csize), csize)
#csize <- csize %>% mutate(imageID = as.numeric(gsub("X", "", imageID)))
#row.names(csize) <- NULL
#rownames(csize) <- csize[,1]
#csize<-csize[c("imageID")]
#largelwing<-read.csv("/Users/Apple/Desktop/largest Lwing.csv", header = F)
largewingrep2d<-largewing2d[rep(seq_len(nrow(largewing2d)), 13283), ]

#largebee2D<-replace(largebee2D,1:13341 ,45.46)
# largebee2D<-replace(bees_L_2D,13342:26682 ,26.94)
# largebee2D<-replace(largebee2D,26683:40023 ,40.184)
 #largebee2D<-replace(largebee2D,40024:53364 ,26.656)
 #largebee2D<-replace(largebee2D,53365:66705 ,40.83)
 #largebee2D<-replace(largebee2D,66706:80046 ,27.65)
 #largebee2D<-replace(largebee2D,80047:93387 ,41.80)
 #largebee2D<-replace(largebee2D,93388:106728 ,29.413)
 #largebee2D<-replace(largebee2D,106729:120069 ,40.49)
 #largebee2D<-replace(largebee2D,120070:133410 ,29.44)
 #largebee2D<-replace(largebee2D,133411:146751 ,39.61)
 #largebee2D<-replace(largebee2D,146752:160092 ,28.89)
 #largebee2D<-replace(largebee2D,160093:173433 ,37.925)
 
rownames<-rownames(bees_L_2D)
#largewingrep2d1<-cbind(largewingrep2d2,rownames)
rownames(largewingrep2d2) <- rownames
#largewingrep2d1 <- largewingrep2d1[,-27]


bees_L_2D <- two.d.array(wingdata)
largewing2d<-two.d.array(largewingdata)
bbind<-rbind(bees_L_2D, largewingrep2d2) %>% arrayspecs(p=13, k=2) 
bbindL<-bees_L_2D %>% arrayspecs(p=13, k=2) 

largewingdata<-readland.tps("/Users/Apple/Desktop/left_2_13_15LM/largewing.txt", specID=c("imageID"))
View(bees_L_2D)

is.numeric(bbind)
is.numeric(lar)
largewingrep2d2<-as.numeric( as.character(largewingrep2d1$V1))
largewingrep2d2<-as.matrix(sapply(largewingrep2d1, as.numeric))  
largewingrep2d2 <- matrix(as.numeric(largewingrep2d1),    # Convert to numeric matrix
                  ncol = ncol(largewingrep2d1))
row1<-rep(rownames,2)
side<-c(rep(1,13283), rep(2,13283))

##procd.lm
datframe<-geomorph.data.frame(bee_gpa, ind= row1 )
fit1<-procD.lm(coords ~ Csize ,
                data = datframe , iter = 999, turbo = TRUE,
                RRPP = F, print.progress = F) # randomize raw values

<-fit1$fitted.values

 rownames(bees_L_2D) = make.names(rownames, unique=TRUE)


bee_gpaL <- gpagen(bbindL,ProcD = T)
plot(bee_gpaL)
pdist <- bee_gpalarge$Csize %>% as.matrix()

csizelarge <- bee_gpalarge$Csize %>% as.matrix()
bee_gpalarge$Csize<-csizelarge
coordslarge <- bee_gpalarge$coords %>% as.matrix()
coordslarge
rownames(coordslarge) = make.names(coordslarge, unique=TRUE)
bee_gpalarge$coords<-coordslarge
plot(bee_gpalarge)

reveal.model.designs(fit1)
anovafit<-anova(fit1, fit2 )
gp <-  interaction(datframe$coords, datframe$Csize)
pw <- pairwise(fit1, groups = datframe$Csize)
plot(fit1)
data(lizards)
land.pairs<-c(2:14,16:28)
dim(land.pairs)<-c(13,2)

gdf <- geomorph.data.frame(shape = bee_gpa$coords,
                           ind = row1)

asym<-bilat.symmetry(datframe$coords,ind = datframe$ind ,object.sym=F, land.pairs = land.pairs, side=side)

csize<-csize%>% distinct(ImageID, .keep_all = TRUE)
datawinged<-merge(csize,rawdata,by="ImageID")
#final dataset=3652 observations
write.csv(file = "/Users/Apple/Desktop/datawinged.csv",datawinged)

rawdata<-read.csv(file = "/Users/Apple/Desktop/Test1Bumblebee_specimen_metadata_ITD_ProcD_Climate_12K.csv", header = T)
rawdata<- rawdata[order(rawdata$"sample_collected_curated"),]
datawinged <- na.omit(datawinged) 

allpr <- prcomp(datawinged[ ,c(2,3,16)], scale = TRUE)
summary(allpr)
biplot(allpr)


datawinged <- cbind(datawinged, allpr$x)

fit1<-lm(ITD~csize, data = datawinged)
summary(fit1)
 cor(datawinged$ITD,datawinged$csize)

 write.csv(file = "/Users/Apple/Desktop/datawinged.csv",datawinged)
 datawinged<-read.csv(file = "/Users/Apple/Desktop/datawinged.csv", header = T)
 
 is.factor(datawinged$yeargroup2)
 datawinged$yeargroup2<-as.factor(datawinged$yeargroup2)
 
 workereshift<- ggplot(subset(datawinged , caste=="worker"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup2)) +
   stat_ellipse(  geom = "polygon", col="black", alpha = 0.5) + ylim(-3,3)+xlim(-4,2)+
   geom_point(shape = 21, col = "black") +ggtitle("Worker ellipse shift over time")

 droneeshift<-ggplot(subset(datawinged , caste=="drone"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup2)) +
   stat_ellipse(  geom = "polygon", col = "black", alpha = 0.5) + ylim(-3,3)+xlim(-4,2)+
   geom_point(shape = 21, col = "black") +ggtitle("Drone ellipse shift over time")

 queeneshift<-ggplot(subset(datawinged , caste=="queen"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup2)) +
   stat_ellipse(  type= "t" ,geom = "polygon", col = "black", alpha = 0.5) + ylim(-3,3)+xlim(-2,4)+
   geom_point(shape = 21, col = "black") +ggtitle("Queen ellipse shift over time") 
 
 ggplot(subset(datawinged , caste=="queen"&sample_species_ID=="lapidarius"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup2)) +
   stat_ellipse( type = "euclid", level =0.5 , geom = "polygon", col = "black", alpha = 0.5) + ylim(-3,3)+xlim(-2,5)+
   geom_point(shape = 21, col = "black") +ggtitle("Queen Hort ellipse shift over time") 

 ##area
 pbw = ggplot_build(workereshift)
 elw = pbw$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well
 
 #Make an empty dataframe where your data is going to end up
 ellipseAreaListw = data.frame()
 
 #Identify unique groups (in my case colour)
 group_list = unique(elw$fill)
 group_list
 #Start the loop
 for (i in 1:length(group_list)) {
   
   #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
   Object = subset(elw, elw$fill == group_list[i])
   Object<-na.omit(Object)
   #Remove grouping column
   Object = Object[-1]
   
   # Center of ellipse
   ctr = MASS::cov.trob(Object)$center  
   
   # Calculate distance to center from each point on the ellipse
   dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
   
   # Calculate area of ellipse from semi-major and semi-minor axes. 
   # These are, respectively, the largest and smallest values of dist2center. 
   value = pi*min(dist2center)*max(dist2center)
   
   #Store in the area list
   ellipseAreaListw = rbind(ellipseAreaListw, data.frame(group_list[i], value))
 }  
 
 pbd = ggplot_build(droneeshift)
 eld = pbd$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well
 
 #Make an empty dataframe where your data is going to end up
 ellipseAreaListd = data.frame()
 
 #Identify unique groups (in my case colour)
 group_list = unique(eld$fill)
 
 #Start the loop
 for (i in 1:length(group_list)) {
   
   #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
   Object = subset(eld, eld$fill == group_list[i])
   Object<-na.omit(Object)
   #Remove grouping column
   Object = Object[-1]
   
   # Center of ellipse
   ctr = MASS::cov.trob(Object)$center  
   
   # Calculate distance to center from each point on the ellipse
   dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
   
   # Calculate area of ellipse from semi-major and semi-minor axes. 
   # These are, respectively, the largest and smallest values of dist2center. 
   value = pi*min(dist2center)*max(dist2center)
   
   #Store in the area list
   ellipseAreaListd = rbind(ellipseAreaListd, data.frame(group_list[i], value))
 }  
 
 
 pbq = ggplot_build(queeneshift)
 elq = pbq$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well
 
 #Make an empty dataframe where your data is going to end up
 ellipseAreaListq = data.frame()
 
 #Identify unique groups (in my case colour)
 group_list = unique(elq$fill)
 group_list
 #Start the loop
 for (i in 1:length(group_list)) {
   
   #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
   Object = subset(elq, elq$fill == group_list[i])
   Object<-na.omit(Object)
   #Remove grouping column
   Object = Object[-1]
   
   # Center of ellipse
   ctr = MASS::cov.trob(Object)$center  
   
   # Calculate distance to center from each point on the ellipse
   dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
   
   # Calculate area of ellipse from semi-major and semi-minor axes. 
   # These are, respectively, the largest and smallest values of dist2center. 
   value = pi*min(dist2center)*max(dist2center)
   
   #Store in the area list
   ellipseAreaListq = rbind(ellipseAreaListq, data.frame(group_list[i], value))
 }  
 
 ellipseAreaListd [ellipseAreaListd$value == 0.8120307 ]<- 2.6269967
 ellipseAreaListd<- ellipseAreaListd [-c(8),]
 ellipseAreaListd$yeargroup<-c(1:7)
ggplot(ellipseAreaListd, aes(yeargroup,value))+geom_line()+ggtitle("Drone ellipse size change")+geom_smooth(method = lm)

ellipseAreaListq<- ellipseAreaListq [-c(8),]
ellipseAreaListq$yeargroup<-c(1:7)
ggplot(ellipseAreaListq, aes(yeargroup,value))+geom_line()+ggtitle("Queen ellipse size change")+ geom_smooth(method = lm)


ellipseAreaList$yeargroup<-c(1:7)
ggplot(ellipseAreaList, aes(yeargroup,value))+geom_line()+ggtitle("Worker ellipse size change")+ geom_smooth(method = lm)
##Time vs spp.
lapeshift<-ggplot(subset(datawinged, sample_species_ID=="lapidarius"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup2)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("Lap ellipse shift over time") 

pblap = ggplot_build(lapeshift)
ellap = pblap$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListlap = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(ellap$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(ellap, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListlap = rbind(ellipseAreaListlap, data.frame(group_list[i], value))
}  
ellipseAreaListlap$yeargroup<-c(1:7)
ggplot(ellipseAreaListlap, aes(yeargroup,value))+geom_line()+ggtitle("Lapidarius ellipse size change")+ geom_smooth(method = lm, level = 0.9)

#mus

#museshift<-
ggplot(subset(datawinged, sample_species_ID=="muscorum"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup2)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("Mus ellipse shift over time") 

pbmus = ggplot_build(museshift)
elmus = pbmus$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListmus = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(elmus$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(elmus, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListmus = rbind(ellipseAreaListmus, data.frame(group_list[i], value))
}  
ellipseAreaListmus$yeargroup<-c(1:6)
ggplot(ellipseAreaListmus, aes(yeargroup,value))+geom_line()+ggtitle("Muscorum ellipse size change")+ geom_smooth(method = lm, level = 0.9)

#Pas
paseshift<-ggplot(subset(datawinged, sample_species_ID=="pascuorum"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup2)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("Pas ellipse shift over time") 

pbpas = ggplot_build(paseshift)
elpas = pbpas$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListpas = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(elpas$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(elpas, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListpas = rbind(ellipseAreaListpas, data.frame(group_list[i], value))
}  
ellipseAreaListpas$yeargroup<-c(1:7)
ggplot(ellipseAreaListpas, aes(yeargroup,value))+geom_line()+ggtitle("Pascuorum ellipse size change")+ geom_smooth(method = lm, level = 0.9)
# syl 
syleshift<-ggplot(subset(datawinged, sample_species_ID=="sylvarum"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup2)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("syl ellipse shift over time") 

pbsyl = ggplot_build(syleshift)
elsyl = pbsyl$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListsyl = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(elsyl$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(elsyl, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListsyl = rbind(ellipseAreaListsyl, data.frame(group_list[i], value))
}  
ellipseAreaListsyl$yeargroup<-c(1:4)
ggplot(ellipseAreaListsyl, aes(yeargroup,value))+geom_line()+ggtitle("Sylvarum ellipse size change")+ geom_smooth(method = lm, level = 0.9)
#hort 


horeshift<-ggplot(subset(datawinged, sample_species_ID=="hortorum"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup2)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("hor ellipse shift over time") 

pbhor = ggplot_build(horeshift)
elhor = pbhor$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListhor = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(elhor$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(elhor, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListhor = rbind(ellipseAreaListhor, data.frame(group_list[i], value))
}  
ellipseAreaListhor$yeargroup<-c(1:7)
ggplot(ellipseAreaListhor, aes(yeargroup,value))+geom_line()+ggtitle("Hortorum ellipse size change")+ geom_smooth(method = lm, level = 0.9)

#for all 
alleshift<-ggplot(datawinged, aes( PC1, PC2, col = yeargroup2, fill = yeargroup2)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("all ellipse shift over time") 

pball = ggplot_build(alleshift)
elall = pball$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListall = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(elall$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(elall, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListall = rbind(ellipseAreaListall, data.frame(group_list[i], value))
}  
ellipseAreaListall$yeargroup<-c(1:7)
ellipseAreaListall$yeargroup<-as.integer(ellipseAreaListall$yeargroup)
ggplot(ellipseAreaListall, aes(yeargroup,value))+geom_line()+ggtitle("All ellipse size change")+ geom_smooth(method = lm, level = 0.9)




#Feb3: re-introduce the 15 year sliding window with 5 year increments
  #create yeargroups in 5 year increments (disused)
fun <- function(x) {
  (x -1895)/15
}

#datawinged$yeargroup3<-c(1)
datawinged$yeargroup10<-lapply(datawinged$sample_collected_curated, fun)
datawinged$yeargroup1910<-as.integer(datawinged$yeargroup10)
 
 ##instead, plot 3 times, 15-year increments but start with different years. 
ellipseAreaListall$yeargroup<-c("1","4","7","10","13","16","19")

#now create shifts starting 1915: line 468~563
fun <- function(x) { (x -1900)/15}

datawinged$yeargroup1915<-lapply(datawinged$sample_collected_curated, fun)
datawinged$yeargroup1915<-as.integer(datawinged$yeargroup1915)
# ggplot by this yeargroup
datawinged$yeargroup1915<-as.character(datawinged$yeargroup1915)
alleshift1915<-ggplot(datawinged, aes( PC1, PC2, col = yeargroup2, fill = yeargroup1915)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("all ellipse shift 1915") 
#obtain ellipse size 

pball1915 = ggplot_build(alleshift1915)
elall1915 = pball1915$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListall1915 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(elall1915$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(elall1915, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListall1915 = rbind(ellipseAreaListall1915, data.frame(group_list[i], value))
}  
ellipseAreaListall1915$yeargroup<-c(NA,2,5,8,11,14,17,NA)
ellipseAreaListall1915<-na.omit(ellipseAreaListall1915)
ellipseAreaListallA<-rbind(ellipseAreaListall1915,ellipseAreaListall)
  #shifts starting 1920
fun <- function(x) { (x -1905)/15}

datawinged$yeargroup1920<-lapply(datawinged$sample_collected_curated, fun)
datawinged$yeargroup1920<-as.integer(datawinged$yeargroup1920)
# ggplot by this yeargroup
datawinged$yeargroup1920<-as.character(datawinged$yeargroup1920)
alleshift1920<-ggplot(datawinged, aes( PC1, PC2, col = yeargroup2, fill = yeargroup1920)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("all ellipse shift 1920") 
#obtain ellipse size 

pball1920 = ggplot_build(alleshift1920)
elall1920 = pball1920$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListall1920 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(elall1920$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(elall1920, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListall1920 = rbind(ellipseAreaListall1920, data.frame(group_list[i], value))
}  

ellipseAreaListall1920$yeargroup<-c(3,6,9,12,15,18)
ellipseAreaListall1920<-na.omit(ellipseAreaListall1920)
ellipseAreaListallA<-rbind(ellipseAreaListall1915,ellipseAreaListall,ellipseAreaListall1920)
  #ggplot the area shifts 
ellipseAreaListallA$yeargroup<-as.integer(ellipseAreaListallA$yeargroup)
ggplot(ellipseAreaListallA, aes(yeargroup,value))+geom_line()+ggtitle("All ellipse size change")+ geom_smooth(method = lm, level = 0.95)

## move on to Sylvarum, the best fit in previous test :line 565~652
ellipseAreaListsyl$yeargroup<-c(1,4,7,10)

fun <- function(x) { (x -1900)/15}

# ggplot by this yeargroup
syleshift1915<-ggplot(subset(datawinged, sample_species_ID=="sylvarum"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup1915)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("syl ellipse shift 1915") 
#obtain ellipse size 

pbsyl1915 = ggplot_build(syleshift1915)
elsyl1915 = pbsyl1915$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListsyl1915 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(elsyl1915$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(elsyl1915, fill == group_list[6])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListsyl1915 = rbind(ellipseAreaListsyl1915, data.frame(group_list[i], value))
}  
ellipseAreaListsyl1915$yeargroup<-c(2)
#ellipseAreaListsylA<-read.csv(file = "/Users/Apple/Desktop/ellipseAreaListsylA.csv")
# ggplot by this yeargroup
syleshift1920<-ggplot(subset(datawinged, sample_species_ID=="sylvarum") ,aes( PC1, PC2, col = yeargroup2, fill = yeargroup1920)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")  
#obtain ellipse size 

pb = ggplot_build(syleshift1920)

el = pb$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListsyl1920 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(el$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(el, fill == group_list[5])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListsyll1920 = rbind(ellipseAreaListsyl1920, data.frame(group_list[i], value))
}  

ellipseAreaListsyll1920$yeargroup<-c(5)
ellipseAreaListsylA<-rbind(ellipseAreaListsyl1915,ellipseAreaListsyl,ellipseAreaListsyll1920)
#ggplot the area shifts 
ellipseAreaListsylA$yeargroup<-as.integer(ellipseAreaListsylA$yeargroup)
ggplot(ellipseAreaListsylA, aes(yeargroup,value))+geom_line()+ggtitle("syl ellipse size change")+ geom_smooth(method = lm, level = 0.95)



##move on to pas: line 652~745
ellipseAreaListpas$yeargroup<-c(1,4,7,10,13,17,20)
view(ellipseAreaListpas)
fun <- function(x) { (x -1900)/15}

# ggplot by this yeargroup
paseshift1915<-ggplot(subset(datawinged, sample_species_ID=="pascuorum"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup1915)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("pas ellipse shift 1915") 
#obtain ellipse size 

pb = ggplot_build(paseshift1915)
el = pb$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListpas1915 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(el$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(el, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListpas1915 = rbind(ellipseAreaListpas1915, data.frame(group_list[i], value))
}  
ellipseAreaListpas1915$yeargroup<-c(NA,2,5,8,11,14,17,NA)
ellipseAreaListpas1915<-na.omit(ellipseAreaListpas1915)
#ellipseAreaListsylA<-read.csv(file = "/Users/Apple/Desktop/ellipseAreaListsylA.csv")
# ggplot by this yeargroup 1920
paseshift1920<-ggplot(subset(datawinged, sample_species_ID=="pascuorum") ,aes( PC1, PC2, col = yeargroup2, fill = yeargroup1920)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")  
#obtain ellipse size 

pb = ggplot_build(paseshift1920)

el = pb$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListpas1920 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(el$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(el, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListpas1920 = rbind(ellipseAreaListpas1920, data.frame(group_list[i], value))
}  

ellipseAreaListpas1920$yeargroup<-c(NA,3,6,9,12,15,18,NA)
ellipseAreaListpas1920<-na.omit(ellipseAreaListpas1920)
ellipseAreaListpasA<-rbind(ellipseAreaListpas1915,ellipseAreaListpas,ellipseAreaListpas1920)
#ggplot the area shifts 
ellipseAreaListpasA$yeargroup<-as.integer(ellipseAreaListpasA$yeargroup)
ggplot(ellipseAreaListpasA, aes(yeargroup,value))+geom_line()+ggtitle("pas ellipse size change")+ geom_smooth(method = lm, level = 0.95)

##move on to mos:line 747~836

ellipseAreaListmus$yeargroup<-c(1,4,7,10,13,16)

# ggplot by this yeargroup
museshift1915<-ggplot(subset(datawinged, sample_species_ID=="muscorum"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup1915)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("mus ellipse shift 1915") 
#obtain ellipse size 

pb = ggplot_build(museshift1915)
el = pb$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListmus1915 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(el$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(el, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListmus1915 = rbind(ellipseAreaListmus1915, data.frame(group_list[i], value))
}  
ellipseAreaListmus1915$yeargroup<-c(NA,2,5,8,11,14,NA)
ellipseAreaListmus1915<-na.omit(ellipseAreaListmus1915)
#ellipseAreaListsylA<-read.csv(file = "/Users/Apple/Desktop/ellipseAreaListsylA.csv")
# ggplot by this yeargroup 1920
##species trait space by pca
museshift1920<-
  ggplot(traitdata2 ,aes( PC1 , ProcD, col = sample_species_ID, fill = sample_species_ID)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) + ylim(-0.03,0.25)+
  geom_point(shape = 21, col = "black")  + theme_bw() + ggtitle("Specis Trait Space by PCA")
  +theme(axis.text=element_text(size=15))
   

colnames(traitdata$sample_species_ID)<- "Species_ID"
traitdata2<-subset(traitdata, traitdata$sample_species_ID !="sylvarum")
#obtain ellipse size 


pb = ggplot_build(museshift1920)

el = pb$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListmus1920 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(el$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(el, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListmus1920 = rbind(ellipseAreaListmus1920, data.frame(group_list[i], value))
}  

ellipseAreaListmus1920$yeargroup<-c(NA,3,6,9,12,15,NA)
ellipseAreaListmus1920<-na.omit(ellipseAreaListmus1920)
ellipseAreaListmusA<-rbind(ellipseAreaListmus1915,ellipseAreaListmus,ellipseAreaListmus1920)
#ggplot the area shifts 
ellipseAreaListmusA$yeargroup<-as.integer(ellipseAreaListmusA$yeargroup)
ggplot(ellipseAreaListmusA, aes(yeargroup,value))+geom_line(size=4)+ggtitle("mus ellipse size change")+ geom_smooth(method = lm, level = 0.95, fill= "blue")+ geom_smooth(fill="orange", level=0.95)


##move on to hor:line 838~925
ellipseAreaListhor$yeargroup<-c(1,4,7,10,13,16,19)

# ggplot by this yeargroup
horeshift1915<-ggplot(subset(datawinged, sample_species_ID=="hortorum"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup1915)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("hor ellipse shift 1915") 
#obtain ellipse size 

pb = ggplot_build(horeshift1915)
el = pb$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListhor1915 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(el$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(el, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListhor1915 = rbind(ellipseAreaListhor1915, data.frame(group_list[i], value))
}  
ellipseAreaListhor1915$yeargroup<-c(NA,2,5,8,11,14,17,NA)
ellipseAreaListhor1915<-na.omit(ellipseAreaListhor1915)
#ellipseAreaListsylA<-read.csv(file = "/Users/Apple/Desktop/ellipseAreaListsylA.csv")
# ggplot by this yeargroup 1920
horeshift1920<-ggplot(subset(datawinged, sample_species_ID=="hortorum") ,aes( PC1, PC2, col = yeargroup2, fill = yeargroup1920)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")  
#obtain ellipse size 

pb = ggplot_build(horeshift1920)

el = pb$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListhor1920 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(el$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(el, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListhor1920 = rbind(ellipseAreaListhor1920, data.frame(group_list[i], value))
}  

ellipseAreaListhor1920$yeargroup<-c(NA,3,6,9,12,15,18,NA)
ellipseAreaListhor1920<-na.omit(ellipseAreaListhor1920)
ellipseAreaListhorA<-rbind(ellipseAreaListhor1915,ellipseAreaListhor,ellipseAreaListhor1920)
#ggplot the area shifts 
ellipseAreaListhorA$yeargroup<-as.integer(ellipseAreaListhorA$yeargroup)
ggplot(ellipseAreaListhorA, aes(yeargroup,value))+geom_line(size=4)+ggtitle("hor ellipse size change")+ geom_smooth(method = lm, level = 0.95, fill= "blue")+ geom_smooth(fill="orange", level=0.95)

##momve on to lap:line 927~1016

ellipseAreaListlap$yeargroup<-c(1,4,7,10,13,16,19)

# ggplot by this yeargroup
lapeshift1915<-ggplot(subset(datawinged, sample_species_ID=="lapidarius"), aes( PC1, PC2, col = yeargroup2, fill = yeargroup1915)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") +ggtitle("lap ellipse shift 1915") 
#obtain ellipse size 

pb = ggplot_build(lapeshift1915)
el = pb$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListlap1915 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(el$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(el, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListlap1915 = rbind(ellipseAreaListlap1915, data.frame(group_list[i], value))
}  
ellipseAreaListlap1915$yeargroup<-c(NA,2,5,8,11,14,17,NA)
ellipseAreaListlap1915<-na.omit(ellipseAreaListlap1915)
#ellipseAreaListsylA<-read.csv(file = "/Users/Apple/Desktop/ellipseAreaListsylA.csv")
# ggplot by this yeargroup 1920
lapeshift1920<-ggplot(subset(datawinged, sample_species_ID=="lapidarius") ,aes( PC1, PC2, col = yeargroup2, fill = yeargroup1920)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")  
#obtain ellipse size 

pb = ggplot_build(lapeshift1920)

el = pb$data[[2]][c("fill", "x","y")] #If you separated by shape as well then include it here as well and then separate by shape as well

#Make an empty dataframe where your data is going to end up
ellipseAreaListlap1920 = data.frame()

#Identify unique groups (in my case colour)
group_list = unique(el$fill)
group_list
#Start the loop
for (i in 1:length(group_list)) {
  
  #Subset to separate into individual ellipse groups - if you used shape in addition to colour you'll want to change this as well
  Object = subset(el, fill == group_list[i])
  Object<-na.omit(Object)
  #Remove grouping column
  Object = Object[-1]
  
  # Center of ellipse
  ctr = MASS::cov.trob(Object)$center  
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  value = pi*min(dist2center)*max(dist2center)
  
  #Store in the area list
  ellipseAreaListlap1920 = rbind(ellipseAreaListlap1920, data.frame(group_list[i], value))
}  

ellipseAreaListlap1920$yeargroup<-c(NA,3,6,9,12,15,18,NA)
ellipseAreaListlap1920<-na.omit(ellipseAreaListlap1920)
ellipseAreaListlapA<-rbind(ellipseAreaListlap1915,ellipseAreaListlap,ellipseAreaListlap1920)
#ggplot the area shifts 
ellipseAreaListlapA$yeargroup<-as.integer(ellipseAreaListlapA$yeargroup)
ggplot(ellipseAreaListlapA, aes(yeargroup,value))+geom_line(size=4)+ggtitle("lap ellipse size change")+ geom_smooth(method = lm, level = 0.95, fill= "blue")+ geom_smooth(fill="orange", level=0.95)

##check data sample size
datasyl<-subset(datawinged, datawinged$sample_species_ID=="sylvarum")
view(datasyl)
summary(datasyl$yeargroup2)
summary(datasyl$yeargroup1915)
summary(datasyl$yeargroup1920)
datamus<-subset(datawinged, datawinged$sample_species_ID=="muscorum")
summary(datamus$yeargroup1920)


##Feb10 1) saturation curve for ellipse size. 
library(Rcpp)
library(hypervolume)
rownames(datawinged)<-datawinged[,1]
pc<-datawinged[c(28:29)]
pchype<-hypervolume(pc)
pchype

  #a) create random selection ,test with mus (line 1033~1088)
pc<-ggplot(datawinged ,aes( PC1, PC2, col = yeargroup2)) +geom_point(shape = 21, col = "black")  
datamus<-subset(datawinged, sample_species_ID=="muscorum")
datamus<-datamus[c(29:30)]
#ranmus <- datamus[sample(1:nrow(datamus), 50,replace=FALSE),] 
##create datasheet by repeatition 

#loop 1: generate x axis--numbers of sample selected, starting with 3 (smallest # to generate circle)
mussat<-data.frame()
for(t in 5:200){
  
ranmus50=data.frame()
#loop 1.1: create datasheet
for(j in 1:50){
  ranmus<-datamus[sample(1:nrow(datamus), t,replace=FALSE),]
  x<-rep(c(j),times=t)
  ranmus1<-cbind (ranmus,x)
  ranmus50<-rbind(ranmus50,ranmus1)
}
# make ellipses based on selection
ranmus50$x<-as.character(ranmus50$x)
ranmusplot<-ggplot(ranmus50, aes(PC1,PC2, fill=x)) +stat_ellipse(geom = "polygon", alpha=0.5)
#obtain ellipse size 
pb = ggplot_build(ranmusplot)
el = pb$data[[1]][c("fill", "x","y")] 
ellipseAreaListranmus = data.frame()
group_list = unique(el$fill)
group_list
#loop 1.2: calculate ellipse size 
for (i in 1:length(group_list)) {
  Object = subset(el, fill == group_list[i])
  Object<-na.omit(Object)
  Object = Object[-1]
    ctr = MASS::cov.trob(Object)$center  
    dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
  value = pi*min(dist2center)*max(dist2center)
    ellipseAreaListranmus = rbind(ellipseAreaListranmus, data.frame(group_list[i], value))
} 
#add up each try and make list
mussat1<-mean(ellipseAreaListranmus$value)
mussat=rbind(mussat,mussat1)
}
mussat$x <- rownames(mussat)
fun <- function(x) {
  (x+4)
}
mussat$x<-as.integer (mussat$x)
mussat$x1<-lapply(mussat$x,fun)
mussat$y<-mussat$X27.0116572483411
plot(mussat$x,mussat$y) 
mussat$x1<-as.integer (mussat$x1)
ggplot(mussat,aes(x1,y))+ geom_point()+stat_smooth(method="lm",formula=  y ~ I(1/x))+ggtitle("Muscorom Ellispe Size Vs Sample Size")+  xlim(5,200)+ xlab("sample size")


#move on with pas
datapas<-subset(datawinged, sample_species_ID=="pascuorum")
datapas<-datapas[c(29:30)]
#ranmus <- datamus[sample(1:nrow(datamus), 50,replace=FALSE),] 
##create datasheet by repeatition 

#loop 1: generate x axis--numbers of sample selected, starting with 3 (smallest # to generate circle)
passat<-data.frame()
for(t in 5:200){
  
  ranpas50=data.frame()
  #loop 1.1: create datasheet
  for(j in 1:100){
    ranpas<-datapas[sample(1:nrow(datapas), t,replace=FALSE),]
    x<-rep(c(j),times=t)
    ranpas1<-cbind (ranpas,x)
    ranpas50<-rbind(ranpas50,ranpas1)
  }
  # make ellipses based on selection
  ranpas50$x<-as.character(ranpas50$x)
  ranpasplot<-ggplot(ranpas50, aes(PC1,PC2, fill=x)) +stat_ellipse(geom = "polygon", alpha=0.5)
  #obtain ellipse size 
  pb = ggplot_build(ranpasplot)
  el = pb$data[[1]][c("fill", "x","y")] 
  ellipseAreaListranpas = data.frame()
  group_list = unique(el$fill)
  group_list
  #loop 1.2: calculate ellipse size 
  for (i in 1:length(group_list)) {
    Object = subset(el, fill == group_list[i])
    Object<-na.omit(Object)
    Object = Object[-1]
    ctr = MASS::cov.trob(Object)$center  
    dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
    value = pi*min(dist2center)*max(dist2center)
    ellipseAreaListranpas = rbind(ellipseAreaListranpas, data.frame(group_list[i], value))
  } 
  #add up each try and make list
  passat1<-mean(ellipseAreaListranpas$value)
  passat=rbind(passat,passat1)
}
passat$x <- rownames(passat)
fun <- function(x) {
  (x+4)
}
passat$x<-as.integer (passat$x)
passat$x1<-lapply(passat$x,fun)
passat$y<-passat$X20.0610305853595
passat$x1<-as.integer (passat$x1)
ggplot(passat,aes(x1,y))+ geom_point()+stat_smooth(method="lm",formula=  y ~ I(1/x))+ggtitle("Pas Ellispe Size Vs Sample Size")+  xlim(5,200)+ xlab("sample size")

#move on with syl

datasyl<-subset(datawinged, sample_species_ID=="sylvarum")
datasyl<-datasyl[c(29:30)]
#ranmus <- datamus[sample(1:nrow(datamus), 50,replace=FALSE),] 
##create datasheet by repeatition 

#loop 1: generate x axis--numbers of sample selected, starting with 3 (smallest # to generate circle)
sylsat<-data.frame()
for(t in 5:100){
  
  ransyl50=data.frame()
  #loop 1.1: create datasheet
  for(j in 1:50){
    ransyl<-datasyl[sample(1:nrow(datasyl), t,replace=F),]
    x<-rep(c(j),times=t)
    ransyl1<-cbind (ransyl,x)
    ransyl50<-rbind(ransyl50,ransyl1)
  }
  # make ellipses based on selection
  ransyl50$x<-as.character(ransyl50$x)
  ransylplot<-ggplot(ransyl50, aes(PC1,PC2, fill=x)) +stat_ellipse(geom = "polygon")
  #obtain ellipse size 
  pb = ggplot_build(ransylplot)
  el = pb$data[[1]][c("fill", "x","y")] 
  ellipseAreaListransyl = data.frame()
  group_list = unique(el$fill)
  group_list
  #loop 1.2: calculate ellipse size 
  for (i in 1:length(group_list)) {
    Object = subset(el, fill == group_list[i])
    Object<-na.omit(Object)
    Object = Object[-1]
    ctr = MASS::cov.trob(Object)$center  
    dist2center <- sqrt(rowSums((t(t(Object)-ctr))^2))
    value = pi*min(dist2center)*max(dist2center)
    ellipseAreaListransyl = rbind(ellipseAreaListransyl, data.frame(group_list[i], value))
  } 
  #add up each try and make list
  sylsat1<-mean(ellipseAreaListransyl$value)
  sylsat=rbind(sylsat,sylsat1)
}
sylsat$x <- rownames(sylsat)
fun <- function(x) {
  (x+4)
}
sylsat$x<-as.integer (sylsat$x)
sylsat$x1<-lapply(sylsat$x,fun)
sylsat$y<-sylsat$X18.482352883628
sylsat$x1<-as.integer (sylsat$x1)
ggplot(sylsat,aes(x1,y))+ geom_point()+stat_smooth(method="lm",formula=  y ~ I(1/x))+ggtitle("Syl Ellispe Size Vs Sample Size")+  xlim(5,100)+ xlab("sample size")

## fixed sample size yeargroups 

##weather data and yeargroup 
  #Temp   
    
dataweather<-datawinged[-c(1),]
dataweather$temp_WIN<-as.numeric(dataweather$temp_WIN)
alltemp$temp<-as.numeric(alltemp$temp)
alltemp$year<-as.numeric(alltemp$year)
ggplot(dataweather,aes(sample_collected_curated,temp_WIN))+ geom_point()+ ylim(-1,8)+ geom_smooth(level=0.95)+ggtitle("Winter Temp VS year")
ggplot(dataweather,aes(sample_collected_curated,temp_SPR))+ geom_point()+ ylim(4,12)+ geom_smooth(level=0.95)+ggtitle("Spring Temp VS year")
  #precipitation
ggplot(dataweather,aes(sample_collected_curated,rain_ANN))+ geom_point()+ geom_smooth(level=0.95)+ggtitle("Rain_ANN VS year")
    #all four seasons temp and ann precipitation 
wintemp<-dataweather[c(19,10)]
x<-rep(c("win"),times=2653)
wintemp<-cbind(x,wintemp)
colnames(wintemp)<-c("x","temp","year")

sprtemp<-dataweather[c(20,10)]
x<-rep(c("spr"),times=2653)
sprtemp<-cbind(x,sprtemp)
colnames(sprtemp)<-c("x","temp","year")

sumtemp<-dataweather[c(21,10)]
x<-rep(c("sum"),times=2653)
sumtemp<-cbind(x,sumtemp)
colnames(sumtemp)<-c("x","temp","year")

auttemp<-dataweather[c(22,10)]
x<-rep(c("aut"),times=2653)
auttemp<-cbind(x,auttemp)
colnames(auttemp)<-c("x","temp","year")

alltemp<-rbind(wintemp,sprtemp,sumtemp,auttemp)
colnames(alltemp)<-c("Season","Temp", "Year")
ggplot(alltemp,aes(Year,Temp,col= Season))+
  geom_point()+ geom_smooth(level=0.95)+ theme_bw() +labs(x = "Year", y = "Temperature")

  #Off-set yeargroup-ness 
alltemp5<-alltemp
fun<-function(x){(x-5)}
alltemp5$year<-as.integer(alltemp5$year)
alltemp5$year<-lapply(alltemp5$year,fun)
alltemp10<-alltemp5
alltemp10$year<-lapply(alltemp10$year,fun)
repalltemp<-rbind(alltemp,alltemp5,alltemp10)
repalltemp$year<-as.integer(repalltemp$year)
ggplot(repalltemp,aes(year,temp,col=x))+geom_point(alpha=0.5,size=1)+ggtitle("Mimic yeargroup Temp")+ geom_smooth()+xlim(1910,2016)
  #yeargroup-precipitate 
rainann<-dataweather[c(28,10)]
rainann5<-rainann
rainann5$sample_collected_curated<-lapply(rainann5$sample_collected_curated,fun)
rainann10<-rainann5
rainann10$sample_collected_curated<-lapply(rainann10$sample_collected_curated,fun)
reprainann<-rbind(rainann,rainann10,rainann5)
reprainann$sample_collected_curated<-as.integer(reprainann$sample_collected_curated)
ggplot(reprainann,aes(sample_collected_curated,rain_ANN))+geom_point(alpha=0.5,size=1)+ggtitle("Mimic yeargroup rain")+ geom_smooth()

    #maybe temp/rain was special for the three spp. 
  #mus
datamus<-subset(datawinged, datawinged$sample_species_ID=="muscorum")

muswin<-datamus[c(19,10)]
x<-rep(c("win"),times=356)
muswin<-cbind(x,muswin)
colnames(muswin)<-c("x","temp","year")
muswin$temp<-as.numeric(muswin$temp)
ggplot(muswin,aes(year,temp))+geom_point(alpha=0.5,size=1)+ggtitle("Mimic yeargroup rain")+ geom_smooth()
    ##off-set
muswin5<-muswin
muswin5$year<-lapply(muswin5$year,fun)
muswin10<-muswin5
muswin10$year<-lapply(muswin10$year,fun)
muswin<-rbind(muswin,muswin5,muswin10)
muswin$year<-as.integer(muswin$year)
ggplot(muswin,aes(year,temp))+geom_point(alpha=0.5,size=1)+ggtitle("Mus rep Winter temp")+ geom_smooth()+xlim(1910,2016)

    ##lap
datalap<-subset(datawinged, datawinged$sample_species_ID=="lapidarius")

lapwin<-datalap[c(19,10)]
x<-rep(c("win"),times=510)
lapwin<-cbind(x,lapwin)
colnames(lapwin)<-c("x","temp","year")
lapwin$temp<-as.numeric(lapwin$temp)
ggplot(lapwin,aes(year,temp))+geom_point(alpha=0.5,size=1)+ggtitle("Mimic yeargroup rain")+ geom_smooth()
##off-set
lapwin5<-lapwin
lapwin5$year<-lapply(lapwin5$year,fun)
lapwin10<-lapwin5
lapwin10$year<-lapply(lapwin10$year,fun)
lapwin<-rbind(lapwin,lapwin5,lapwin10)
lapwin$year<-as.integer(lapwin$year)
ggplot(lapwin,aes(year,temp))+geom_point(alpha=0.5,size=1)+ggtitle("lap rep Winter temp")+ geom_smooth()+xlim(1910,2016)

    ##hor
datahor<-subset(datawinged, datawinged$sample_species_ID=="hortorum")

horwin<-datahor[c(19,10)]
x<-rep(c("win"),times=573)
horwin<-cbind(x,horwin)
colnames(horwin)<-c("x","temp","year")
horwin$temp<-as.numeric(horwin$temp)
ggplot(horwin,aes(year,temp))+geom_point(alpha=0.5,size=1)+ggtitle("Mimic yeargroup rain")+ geom_smooth()
##off-set
horwin5<-horwin
horwin5$year<-lapply(horwin5$year,fun)
horwin10<-horwin5
horwin10$year<-lapply(horwin10$year,fun)
horwin<-rbind(horwin,horwin5,horwin10)
horwin$year<-as.integer(horwin$year)
ggplot(horwin,aes(year,temp))+geom_point(alpha=0.5,size=1)+ggtitle("hor rep Winter temp")+ geom_smooth()+xlim(1910,2016)

    #precipitation
horrain<-datahor[c(28,10)]
colnames(horrain)<-c("rain","year")
ggplot(horrain,aes(year,rain))+geom_point(alpha=0.5,size=1)+ggtitle("hor Mimic yeargroup rain")+ geom_smooth()
##off-set
horrain5<-horrain
horrain5$year<-lapply(horrain5$year,fun)
horrain10<-horrain5
horrain10$year<-lapply(horrain10$year,fun)
horrain<-rbind(horrain,horrain5,horrain10)
horrain$year<-as.integer(horrain$year)
ggplot(horrain,aes(year,rain))+geom_point(alpha=0.5,size=1)+ggtitle("hor rep rain")+ geom_smooth()+xlim(1910,2016)

  #lap
laprain<-datalap[c(28,10)]
colnames(laprain)<-c("rain","year")
ggplot(laprain,aes(year,rain))+geom_point(alpha=0.5,size=1)+ggtitle("hor Mimic yeargroup rain")+ geom_smooth()
##off-set
laprain5<-laprain
laprain5$year<-lapply(laprain5$year,fun)
laprain10<-laprain5
laprain10$year<-lapply(laprain10$year,fun)
laprain<-rbind(laprain,laprain5,laprain10)
laprain$year<-as.integer(laprain$year)
ggplot(laprain,aes(year,rain))+geom_point(alpha=0.5,size=1)+ggtitle("lap rep rain")+ geom_smooth()+xlim(1910,2000)

  #mus
musrain<-datamus[c(28,10)]
colnames(musrain)<-c("rain","year")
ggplot(musrain,aes(year,rain))+geom_point(alpha=0.5,size=1)+ggtitle("hor Mimic yeargroup rain")+ geom_smooth()
##off-set
musrain5<-musrain
musrain5$year<-lapply(musrain5$year,fun)
musrain10<-musrain5
musrain10$year<-lapply(musrain10$year,fun)
musrain<-rbind(musrain,musrain5,musrain10)
musrain$year<-as.integer(musrain$year)
ggplot(musrain,aes(year,rain))+geom_point(alpha=0.5,size=1)+ggtitle("mus rep rain")+ geom_smooth()+xlim(1910,2000)


#3. Compare pre1960 and post 1960 mus
datawinged$yeargroup1960 <- ifelse(datawinged$sample_collected_curated >=1960, "2", "1")
dataqeen<-subset(datawinged, datawinged$caste=="queen")
datadrone<-subset(datawinged ,datawinged$caste=="drone")
datawork<-subset(datawinged, datawinged$caste=="worker")
ggplot(datawork, aes( PC1, PC2, col = yeargroup1960, fill = yeargroup1960)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") 
datahorq<-subset(dataqeen,dataqeen$sample_species_ID=="hortorum")
datahorw<-subset(datawork,dataqeen$sample_species_ID=="hortorum")
datahord<-subset(datadrone,dataqeen$sample_species_ID=="hortorum")
ggplot(datahorq, aes( PC1, PC2, col = yeargroup1960, fill = yeargroup1960)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") + ggtitle("hor queen ")

datawinged$trigroup<-ifelse(datawinged$sample_collected_curated<= 1940,"1",
                             ifelse(datawinged$sample_collected_curated>1940 & datawinged$sample_collected_curated<=1970,"2",
                                    ifelse(datawinged$sample_collected_curated>1970,"3","4")))

ggplot(datawinged, aes( PC1, PC2, fill = trigroup)) +
  stat_ellipse( type = "t", geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black") 


##Feb 17 
#Task1 investigate how PCA ellipses are generated, and its power

library(BiodiversityR)

traitdata<-datawinged[c(1:3,5,7,9,10,17)]
rownames(traitdata)<-traitdata[,1]
traitdata<-traitdata[-c(1)]
PCAsignificance(pca1)
ordiplot(allpr, type = "text")
library(ggvegan)
library(vegan)
 autoplot(pca1)
datapca<-traitdata[,c(2,1,6,)]
pca1<-rda(datapca, scale = T)
summary(pca1)
plot(pca1, xlim=c(-1,1),ylim=c(-1.5,0.45))  
pcplot1<-ordiplot(pca1, type="text")
ordiellipse(pcplot1, groups = traitdata$sample_species_ID)
ordiellipse(pca1, traitdata$yeargroup2 ,display = "sites", kind="sd", conf = 0.9, col = traitdata$yeargroup2)

##hypervolume 
library(hypervolume)
hypedata<-traitdata[,c(3,11)]
hv<-hypervolume(hypedata)
summary(hv)
plot(hv)

partpr<-prcomp(traitdata[ ,c(2,8)], scale = TRUE)
traitdata<-cbind(traitdata,partpr$x)
summary(partpr)
## Mar3 
library(alphahull)
##1. create hypervolume for pas, yeargroup 1 only 
pasyg1<-subset(traitdata, yeargroup2=="1" & sample_species_ID  =="pascuorum")

pasyg1<-pasyg1[,c(3,11)]
ranpas1<-pasyg1[sample(1:nrow(pasyg1), 100,replace=FALSE),]

pasyg1hv = hypervolume_svm(ranpas1, name = 'pas1')
get_centroid_weighted(pasyg1hv)
get_volume(pasyg1hv)
## check pas yg=6 
pasyg6<-subset(traitdata, yeargroup2=="6" & sample_species_ID  =="pascuorum")
pasyg6<-pasyg6[,c(3,11)]
ranpas6<-pasyg6[sample(1:nrow(pasyg6), 100,replace=T),]

pasyg6hv = hypervolume_svm(ranpas6, name = 'pas6')
get_centroid_weighted(pasyg1hv)
get_volume(pasyg6hv)
plot(pasyg1hv)


paslist<-hypervolume_join(pasyg1hv,pasyg6hv)
plot.HypervolumeList(paslist, contour.type="kde",show.contour = F)

## check all spp yg1 vs yg 6
allyg1<-subset(traitdata, yeargroup2=="1"| yeargroup2=="2")
allyg1<-allyg1[,c(3,11)]
allyg1hv = hypervolume_svm(allyg1, name = 'all1')

allyg6<-subset(traitdata, yeargroup2== "5"| yeargroup2=="6"| yeargroup2=="7")
allyg6<-allyg6[,c(3,11)]
allyg6hv = hypervolume_svm(allyg6, name = 'all6')
alllist<-hypervolume_join(allyg1hv,allyg6hv)
plot(alllist)

#test for min sample size 
test.size<-data.frame()
for(t in 5:200){
  
  test50=data.frame()
  #loop 1.1: create datasheet
  for(j in 1:50){
    ranmus<-hypedata[sample(1:nrow(datamus), t,replace=FALSE),]
    x<-rep(c(j),times=t)
    ranmus1<-cbind (ranmus,x)
    test50<-rbind(test50,ranmus1)
  }
  # make ellipses based on selection
  test50<-test50[,c(1:2)]
  hv50<-hypervolume_svm(test50)
  plot(hv50)
  vol50<-get_volume(hv50)
    #obtain ellipse size 




  vollist = rbind(vollist, vol50)
  
  #add up each try and make list
  mussat1<-mean(vollist$X0.939320100823515)
  test.size=rbind(test.size,mussat1)
}
test.size$x <- rownames(test.size)
fun <- function(x) {
  (x+4)
}
test.size$x<-as.integer (test.size$x)

test.size$x<-as.integer (test.size$x)
ggplot(test.size,aes(x,X0.939320100823515))+ geom_point()+stat_smooth(method="lm",formula=  y ~ I(1/x))+ggtitle("volume vs sample size")+  
  xlim(5,200)+ xlab("Sample Size") + ylab("Volume")+ theme_bw()

##move on to pas 7 yeargroups 
datawinged<-cbind(datawinged,partpr$x)
traitdata<-datawinged[,c(1,3,5,7,10,33,34,35)]
datapas<-subset(traitdata, sample_species_ID=="pascuorum")

##first generate 7 yeargroups based on "yeargroup2" 

pascent<-data.frame()
pascentmean<-data.frame()
pascentmean2<-data.frame()
for (i in 1:7){
  pasdata<-subset(datapas, yeargroup2== i )
  pasdata<-pasdata[,c(8,2)]
  for (j in 1:10){
    hv = hypervolume_svm(pasdata)
    cent<-get_centroid(hv)
    pascent<-rbind(pascent,cent)
    colnames(pascent)<-c("x","y")
  }    
  pascentmean<-cbind(mean(pascent$x), mean(pascent$y))
  pascent<-data.frame()
  pascentmean2<-rbind(pascentmean2,pascentmean)
  pascentmean<-data.frame()
}
  
pascentmean2$x<-rownames(pascentmean2)
ggplot(pascentmean2, aes(V1,V2)) + geom_point()+geom_path()+ggtitle("Pas centroid shift")

 ##muscorum?
datamus<-subset(traitdata, sample_species_ID=="muscorum")

muscent<-data.frame()
muscentmean<-data.frame()
muscentmean2<-data.frame()
for (i in 1:5){
  musdata<-subset(datamus, yeargroup2== i )
  musdata<-musdata[,c(2,8)]
  for (j in 1:50){
    musdata2<-musdata[sample(1:nrow(musdata), 100,replace=T),]
    hv = hypervolume_svm(musdata2)
    cent<-get_centroid(hv)
    muscent<-rbind(muscent,cent)
    colnames(muscent)<-c("x","y")
    muscentmean<-cbind(mean(muscent$x), mean(muscent$y))
    muscent<-data.frame()
  }
  muscentmean2<-rbind(muscentmean2,muscentmean)
}

muscentmean2$x<-rownames(muscentmean2)
ggplot(muscentmean2, aes(V1,V2)) + geom_point()+geom_path()+ggtitle("Mus centroid shift")


##hort
datahor<-subset(traitdata, sample_species_ID=="hortorum")

muscent<-data.frame()
muscentmean<-data.frame()
muscentmean2<-data.frame()
for (i in 1:6){
  musdata<-subset(datahor, yeargroup2== i )
  musdata<-musdata[,c(2,8)]
  for (j in 1:50){
    musdata2<-musdata[sample(1:nrow(musdata), 100,replace=T),]
    hv = hypervolume_svm(musdata2)
    cent<-get_centroid(hv)
    muscent<-rbind(muscent,cent)
    colnames(muscent)<-c("x","y")
    muscentmean<-cbind(mean(muscent$x), mean(muscent$y))
    muscent<-data.frame()
  }
  muscentmean2<-rbind(muscentmean2,muscentmean)
}

muscentmean2$x<-rownames(muscentmean2)
ggplot(muscentmean2, aes(V1,V2)) + geom_point()+geom_path()+ggtitle("Hor centroid shift")

##lap
datalap<-subset(traitdata, sample_species_ID=="lapidarius")

muscent<-data.frame()
muscentmean<-data.frame()
muscentmean2<-data.frame()
for (i in 1:6){
  musdata<-subset(datalap, yeargroup2== i )
  musdata<-musdata[,c(2,8)]
  for (j in 1:50){
    musdata2<-musdata[sample(1:nrow(musdata), 100,replace=T),]
    hv = hypervolume_svm(musdata2)
    cent<-get_centroid(hv)
    muscent<-rbind(muscent,cent)
    colnames(muscent)<-c("x","y")
    muscentmean<-cbind(mean(muscent$x), mean(muscent$y))
    muscent<-data.frame()
  }
  muscentmean2<-rbind(muscentmean2,muscentmean)
}

muscentmean2$x<-rownames(muscentmean2)
ggplot(muscentmean2, aes(V1,V2)) + geom_point()+geom_path()+ggtitle("Lap centroid shift")


##volume/area shift in lap 
datalap<-subset(traitdata, sample_species_ID=="lapidarius")

hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:6){
  hvdata<-subset(datalap, yeargroup2== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:50){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvolmean2$x<-rownames(hvvolmean2)
hvvolmean2$x<-as.numeric(hvvolmean2$x)
ggplot(hvvolmean2, aes(i,y)) + geom_point()+geom_line()+ggtitle("Lap area shift")

##hor area shift

hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:6){
  hvdata<-subset(datahor, yeargroup2== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:50){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvolmean2$x<-rownames(hvvolmean2)
hvvolmean2$x<-as.numeric(hvvolmean2$x)
ggplot(hvvolmean2, aes(i,y)) + geom_point()+geom_line()+ggtitle("hor area shift")+ theme_bw()
##pas 

hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:6){
  hvdata<-subset(datapas, yeargroup2== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:50){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvolmean2$x<-rownames(hvvolmean2)
hvvolmean2$x<-as.numeric(hvvolmean2$x)
ggplot(hvvolmean2, aes(i,y)) + geom_point()+geom_line()+ggtitle("pas area shift")

##mus
hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:6){
  hvdata<-subset(datamus, yeargroup2== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:50){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

#hvvolmean2$x<-rownames(hvvolmean2)
#hvvolmean2$x<-as.numeric(hvvolmean2$x)
ggplot(hvvolmean2, aes(i,y)) + geom_point()+geom_line()+ggtitle("mus area shift")+theme_bw()+ xlab("Year-group")+ylab("Volume")


##sliding window: mus. By creating 3 * areashifts and assign them yeargroups. 

hvvolmean2["i"][hvvolmean2["i"] == 4] <- 10
hvvolmean2["i"][hvvolmean2["i"] == 3] <- 7
hvvolmean2["i"][hvvolmean2["i"] == 2] <- 4
hvvolmean2["i"][hvvolmean2["i"] == 5] <- 13
hvvolmean2["i"][hvvolmean2["i"] == 6] <- 16
hvvol1<-hvvolmean2

hvvol2<-hvvolmean2
hvvol2["i"][hvvol2["i"] == 6] <- 17
hvvol2["i"][hvvol2["i"] == 5] <- 14
hvvol2["i"][hvvol2["i"] == 4] <- 11
hvvol2["i"][hvvol2["i"] == 3] <- 8
hvvol2["i"][hvvol2["i"] == 2] <- 5
hvvol2["i"][hvvol2["i"] == 1] <- 2

hvvol3<-hvvolmean2
hvvol3["i"][hvvol3["i"] == 6] <- 18
hvvol3["i"][hvvol3["i"] == 5] <- 15
hvvol3["i"][hvvol3["i"] == 4] <- 12
hvvol3["i"][hvvol3["i"] == 3] <- 9
hvvol3["i"][hvvol3["i"] == 2] <- 6
hvvol3["i"][hvvol3["i"] == 1] <- 3

hvvolA<-rbind(hvvol1,hvvol2,hvvol3)

ggplot(hvvolA, aes(i,y)) + geom_point()+geom_line()+ggtitle("Muscorum Volume Over Time")+theme_bw()+ xlab("Year-group")+ylab("Volume")
write.csv(file = "/Users/Apple/Desktop/hvvolA mus.csv",hvvolA)
##move on to pas
#loop1 creates yg 1,4,7,10,13,16
hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:6){
  hvdata<-subset(datapas, yeargroup2== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:10){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvolmean2["i"][hvvolmean2["i"] == 4] <- 10
hvvolmean2["i"][hvvolmean2["i"] == 3] <- 7
hvvolmean2["i"][hvvolmean2["i"] == 2] <- 4
hvvolmean2["i"][hvvolmean2["i"] == 5] <- 13
hvvolmean2["i"][hvvolmean2["i"] == 6] <- 16
hvvol1<-hvvolmean2

#loop2 creates yg 2,5,8,11,14,17
hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:5){
  hvdata<-subset(traitdata, yeargroup1915== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:10){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvol2<-hvvolmean2
hvvol2["i"][hvvol2["i"] == 5] <- 14
hvvol2["i"][hvvol2["i"] == 4] <- 11
hvvol2["i"][hvvol2["i"] == 3] <- 8
hvvol2["i"][hvvol2["i"] == 1] <- 5
hvvol2["i"][hvvol2["i"] == 1] <- 2

#loop3 creates 3,6,9,12,15
hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:5){
  hvdata<-subset(traitdata, yeargroup1920== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:10){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvol3<-hvvolmean2
hvvol3["i"][hvvol3["i"] == 5] <- 15
hvvol3["i"][hvvol3["i"] == 4] <- 12
hvvol3["i"][hvvol3["i"] == 3] <- 9
hvvol3["i"][hvvol3["i"] == 1] <- 6
hvvol3["i"][hvvol3["i"] == 1] <- 3

hvvolA<-rbind(hvvol1,hvvol2,hvvol3)

ggplot(hvvolA, aes(i,y)) + geom_point()+geom_line()+ggtitle("Pas Volume Over Time")+ xlab("Year-group")+ylab("Volume")
write.csv(file = "/Users/Apple/Desktop/hvvolA1 all.csv",hvvolA)

## lap
#loop1 creates yg 1,4,7,10,13,16
hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:6){
  hvdata<-subset(datalap, yeargroup2== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:50){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvolmean2["i"][hvvolmean2["i"] == 4] <- 10
hvvolmean2["i"][hvvolmean2["i"] == 3] <- 7
hvvolmean2["i"][hvvolmean2["i"] == 2] <- 4
hvvolmean2["i"][hvvolmean2["i"] == 5] <- 13
hvvolmean2["i"][hvvolmean2["i"] == 6] <- 16
hvvol1<-hvvolmean2

#loop2 creates yg 2,5,8,11,14,17
hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:5){
  hvdata<-subset(datalap, yeargroup1915== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:50){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvol2<-hvvolmean2
hvvol2["i"][hvvol2["i"] == 5] <- 14
hvvol2["i"][hvvol2["i"] == 4] <- 11
hvvol2["i"][hvvol2["i"] == 3] <- 8
hvvol2["i"][hvvol2["i"] == 2] <- 5
hvvol2["i"][hvvol2["i"] == 1] <- 2

#loop3 creates 3,6,9,12,15
hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:5){
  hvdata<-subset(datalap, yeargroup1920== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:50){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvol3<-hvvolmean2
hvvol3["i"][hvvol3["i"] == 5] <- 15
hvvol3["i"][hvvol3["i"] == 4] <- 12
hvvol3["i"][hvvol3["i"] == 3] <- 9
hvvol3["i"][hvvol3["i"] == 2] <- 6
hvvol3["i"][hvvol3["i"] == 1] <- 3

hvvolA<-rbind(hvvol1,hvvol2,hvvol3)

ggplot(hvvolA, aes(i,y)) + geom_point()+geom_line()+ggtitle("lap area shift")+theme_bw()+ xlab("Year-group")+ylab("Volume")

write.csv(file = "/Users/Apple/Desktop/hvvolA1 lap.csv",hvvolA)

##hor

#loop1 creates yg 1,4,7,10,13,16
hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:6){
  hvdata<-subset(datahor, yeargroup2== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:50){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvolmean2["i"][hvvolmean2["i"] == 4] <- 10
hvvolmean2["i"][hvvolmean2["i"] == 3] <- 7
hvvolmean2["i"][hvvolmean2["i"] == 2] <- 4
hvvolmean2["i"][hvvolmean2["i"] == 5] <- 13
hvvolmean2["i"][hvvolmean2["i"] == 6] <- 16
hvvol1<-hvvolmean2

#loop2 creates yg 2,5,8,11,14,17
hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:5){
  hvdata<-subset(datahor, yeargroup1915== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:50){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvol2<-hvvolmean2
hvvol2["i"][hvvol2["i"] == 5] <- 14
hvvol2["i"][hvvol2["i"] == 4] <- 11
hvvol2["i"][hvvol2["i"] == 3] <- 8
hvvol2["i"][hvvol2["i"] == 2] <- 5
hvvol2["i"][hvvol2["i"] == 1] <- 2

#loop3 creates 3,6,9,12,15
hvvol<-data.frame()
hvvolmean<-data.frame()
hvvolmean2<-data.frame()
for (i in 1:5){
  hvdata<-subset(datahor, yeargroup1920== i )
  hvdata<-hvdata[,c(2,8)]
  for (j in 1:50){
    hvdata2<-hvdata[sample(1:nrow(hvdata), 100,replace=T),]
    hv = hypervolume_svm(hvdata)
    vol<-get_volume(hv)
    hvvol<-rbind(hvvol,vol)
    colnames(hvvol)<-c("y")
    hvvol<-cbind(hvvol,i)
    hvvolmean<-rbind(hvvolmean,hvvol)
    hvvol<-data.frame()
  }
  hvvolmean2<-rbind(hvvolmean2,hvvolmean)
}

hvvol3<-hvvolmean2
hvvol3["i"][hvvol3["i"] == 5] <- 15
hvvol3["i"][hvvol3["i"] == 4] <- 12
hvvol3["i"][hvvol3["i"] == 3] <- 9
hvvol3["i"][hvvol3["i"] == 2] <- 6
hvvol3["i"][hvvol3["i"] == 1] <- 3

hvvolA<-rbind(hvvol1,hvvol2,hvvol3)

ggplot(hvvolA, aes(i,y)) + geom_point()+geom_line()+ggtitle("hor area shift")+ theme_bw()+ xlab("Year-group")+ylab("Volume")
write.csv(file = "/Users/Apple/Desktop/hvvolA1 hor.csv",hvvolA)

##mar 10

# plot 4 spp. in one scatterplot/hypervolume 


allyg5<-subset(traitdata, yeargroup2== "5")
allyg5<-allyg5[,c(2,8)]
allyg5hv = hypervolume_svm(allyg5, name = 'all5')

#create hv list for all 4 spp. 
pasyg5<-subset(allyg5, sample_species_ID == "pascuorum")
pasyg5<-pasyg5[,c(2,8)]
pasyg5hv = hypervolume_svm(pasyg5, name = 'pas5')

musyg5<-subset(allyg5, sample_species_ID == "muscorum")
musyg5<-musyg5[,c(2,8)]
musyg5hv = hypervolume_svm(musyg5, name = 'mus5')

lapyg5<-subset(allyg5, sample_species_ID == "lapidarius")
lapyg5<-lapyg5[,c(2,8)]
lapyg5hv = hypervolume_svm(lapyg5, name = 'lap5')

horyg5<-subset(allyg5, sample_species_ID == "hortorum")
horyg5<-horyg5[,c(2,8)]
horyg5hv = hypervolume_svm(horyg5, name = 'hor5')

alllist5<-hypervolume_join(allyg5hv,pasyg5hv,musyg5hv,lapyg5hv,horyg5hv)
plot.HypervolumeList(alllist5, contour.type="kde",show.contour = T)
get_volume(alllist5)

for ( i in 6){
  allygi<-subset(traitdata, yeargroup2== i)
  allygi<-allygi[,c(2,8)]
  allygihv = hypervolume_svm(allygi, name = 'alli')
  allygi<-subset(traitdata, yeargroup2== i)
  
  #create hv list for all 4 spp. 
  pasygi<-subset(allygi, sample_species_ID == "pascuorum")
  pasygi<-pasygi[,c(2,8)]
  pasygihv = hypervolume_svm(pasygi, name = 'pasi')
  
  musygi<-subset(allygi, sample_species_ID == "muscorum")
  musygi<-musygi[,c(2,8)]
  musygihv = hypervolume_svm(musygi, name = 'musi')
  
  lapygi<-subset(allygi, sample_species_ID == "lapidarius")
  lapygi<-lapygi[,c(2,8)]
  lapygihv = hypervolume_svm(lapygi, name = 'lapi')
  
  horygi<-subset(allygi, sample_species_ID == "hortorum")
  horygi<-horygi[,c(2,8)]
  horygihv = hypervolume_svm(horygi, name = 'hori')
  
  alllisti<-hypervolume_join(allygihv,pasygihv,musygihv,lapygihv,horygihv)
  plot.HypervolumeList(alllisti, contour.type="kde",show.contour = T)
  
}
  
hvlist<-hypervolume_set(allygihv,pasygihv, check.memory = F)
vol6<-get_volume(alllisti)
#volume<- data.frame()
volume <- cbind(volume,vol6)
ggplot(volume, aes(volume)) + geom_point()
is.data.frame(volume)
names(volume)[3] <- "vol3"
write.csv(file = "/Users/Apple/Desktop/volume.csv",volume)
cor(datawinged$csize,datawinged$ITD)
hypervolume_overlap_statistics(hvlist)
cent<-data.frame()
plol<-data.frame()
for ( i in 1:6){
  allygi<-subset(traitdata, yeargroup2== 1)
  #allygi<-traitdata
  
  allygi<-allygi[,c(2,8)]
  allygihv = hypervolume_svm(allygi, name = 'alli')
  allygi<-subset(traitdata, yeargroup2== 1)
  
  #create hv list for all 4 spp. 
  pasygi<-subset(allygi, sample_species_ID == "pascuorum")
  pasygi<-pasygi[,c(2,8)]
  pasygihv = hypervolume_svm(pasygi, name = 'pascuorum')
  
  musygi<-subset(allygi, sample_species_ID == "muscorum")
  musygi<-musygi[,c(2,8)]
  musygihv = hypervolume_svm(musygi, name = 'muscorum')
  
  lapygi<-subset(allygi, sample_species_ID == "lapidarius")
  lapygi<-lapygi[,c(2,8)]
  lapygihv = hypervolume_svm(lapygi, name = 'lapidarius')
  
  horygi<-subset(allygi, sample_species_ID == "hortorum")
  horygi<-horygi[,c(2,8)]
  horygihv = hypervolume_svm(horygi, name = 'hortorum')
  
  alllisti<-hypervolume_join(pasygihv,musygihv,lapygihv,horygihv)
  plot(alllisti, show.contour = F, contour.type = "kde",contour.kde.level = 0.0001 , colors = l)
  cent0<-get_centroid(alllisti)
  cent<-list(cent, cent0)
  hvlist<-hypervolume_set(musygihv, horygihv, check.memory = F)
  paslapol<-hypervolume_overlap_statistics(hvlist)
  plol<-rbind(plol,paslapol)
}
plot(cent)

y<-c(-4,4)
x<-c(0,0.4)
t<-list(x,y)

l<-c(col,col2,col3,col4)
col=colors()[31]
col2=colors()[10]
col3=colors()[87]
col4=colors()[57]


plot(alllisti, show.contour = T, colors = l, limit  = t)



ploverlap<-data.frame()
ploverlap<-plol$X0.489935462131052 #hor

ploverlap<-rbind(plol$X0.462653348777425, ploverlap) #mus

ploverlap<-rbind(plol$X0.506541308365839, ploverlap) #lap

datawinged<-read.csv(file = "/Users/Apple/Desktop/datawinged.csv")
row.names(ploverlap)<-ploverlap$V1
ploverlap<-ploverlap[,-1]
##maps 
datapas1<-subset(datawinged, datawinged$sample_species_ID == "pascuorum"&datawinged$yeargroup2 == 1)
pasyg1sites<-data.frame(longitude = datapas1$long, latitude = datapas1$lat)

ggplot(data = world) +
  geom_sf() +
  geom_point(data = pasyg1sites, aes(x = longitude, y = latitude), size = 1, 
             shape = 23, fill = "darkred") +
  coord_sf(xlim = c(-10,3), ylim = c(50, 63), expand = T)
library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

datapas4<-subset(datawinged, datawinged$sample_species_ID == "hortorum"&datawinged$yeargroup2 == 5)
pasyg4sites<-data.frame(longitude = datapas4$long, latitude = datapas4$lat)
ggplot(data = world) +
  geom_sf() +
  geom_point(data = pasyg4sites, aes(x = longitude, y = latitude), size = 1, 
             shape = 23, fill = "darkred") +
  coord_sf(xlim = c(-10,2), ylim = c(50, 61.5), expand = T)+ ggtitle("Year-group E Hortorum Sample Site")


