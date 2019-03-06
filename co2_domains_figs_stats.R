#### R SCRIPT TO REPRUDUCE THE RESULTS IN THE PUBLICATION: 
#### LANDSCAPE PROCESS DOMAINS DRIVE PATTERNS OF CO2 EVASION FROM RIVER NETWORKS
#### In review in Limnology and Oceanography Letters
#### AUTHORS:  Gerard Rocher-Ros, Ryan Sponseller, William Lidber2, Carl-Magnus Mörth, Reiner Giesler
#### Author of the script and contact email: Gerard Rocher-Ros (g.rocher.ros@gmail.com)
####
#### Last edit: 2019-03-06
#### Before running the script, download the datasets, deposited in the Swedish National Data Service:
#### Dataset Miellajokka can be found here: https://doi.org/10.5878/pxa2-vy55 
#### Global dataset can be found here: https://doi.org/10.5878/77ps-4f21
####
#### The final publication figures and pannels were prepared using Adobe Illustrator
####################################################################################

#### Load  and install packages ####
#List of all packages needed
package_list <- c('ggplot2', 'dplyr', 'readr', 'ggpubr', 'fitdistrplus', 'scales', 'rworldmap',
                  'reshape2','lmSupport')

# Check if there are any packacges missing
packages_missing <- setdiff(package_list, rownames(installed.packages()))

#If we find a package missing, install them
if(length(packages_missing) >= 1) install.packages(packages_missing) 

#Now load all the packages
lapply(package_list, require, character.only = TRUE)

#Set the working directory, you need to change this for he folder where yo have the files
setwd("C:/Users/User/Downloads")

#### Global functions and ggplot theme set ####
## Function to calculate kco2, from Raymond et al., 2012
kCO2 <- function(k600, temp){
  kCO2= k600/ (600/( 1911.1 - 118.11*temp + 3.4527*temp^2 - 0.04132*temp^3 ))^(0.5)
}

# General theme for ggplot
theme_set(theme_classic())

## ggplot parameter to make a linear model, add R2 and add equation
gglm = list(
  geom_smooth(method="lm",color="gray8"),
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.y.npc = 0.9), #put R2 and label
  stat_regline_equation( label.y.npc = 1))


#### Figure 1: Global dataset ####
#### The global dataset and metadata can be found in this repository: https://doi.org/10.5878/77ps-4f21 
#Read file
co2_glob <- read_csv("global_pco2_k.csv") 


### Figure 1a,  export as 6x4 inch and pdf
ggplot(data=co2_glob)+
  geom_point(aes(x=k600,y=pco2), alpha=.5, size=3 )+
  geom_hline(yintercept=400, linetype=2)+
  scale_x_continuous(limits = c(0, 120)) +
  labs(x=expression(k[600]~(m~d^-1)),y=expression(pCO[2]~(ppm)))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14)) 


#Figure 1b: make a global map
w_map <- fortify(map_data('world'))

ggplot()+
  geom_map(data=w_map, map=w_map,
           aes(x=long, y=lat, group=group, map_id=region),
           fill="white", colour="gray40", size=0.5)+
  geom_point(data=co2_glob, aes(Longitude,Latitude), size=1, alpha=.6 , shape=21, fill="gray20", colour="black")+
  geom_point(aes(18.25,68.27), col="red")+
  scale_y_continuous(limits=c(-65,80))+
  labs(y="Latitude", x="Longitude")+
  theme(axis.text = element_text(size=13), axis.title = element_text(size=13))
##export as 3x5 inch pdf


#Figure 1c and MLR. We need to log transform and standardize the dataset
co2glob_log <- co2_glob %>% transmute(STAT_ID, Latitude, Longitude, fco2= log(fco2_gm2, base=10), 
                                       k600=log(kco2, base=10), pco2=log(pco2-380, base=10)) %>% na.omit() %>% 
  filter(pco2>-1000, fco2>-1000) %>% #After the log transform there are many values -Inf, I remove them
  mutate(fco2s=base::scale(fco2), pco2s=base::scale(pco2), k600s=base::scale(k600))


#Figure 1c, export as 6x4 inch and pdf
ggplot(data=co2glob_log, aes(pco2, fco2))+
  geom_point(aes( color=k600), alpha=.5, size=4)+
  scale_color_viridis_c()+
  gglm+
  labs(x=expression(log[10](pCO[2]~(ppm))), y=expression(log[10](F[CO2]~(gC~m^-2~d^-1))))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14)) 

#perform a MLR with the standardized variables
mlr_glob <- lm(fco2s~pco2s+k600s, data=co2glob_log)

#coefficients of the MLR
summary(mlr_glob)

#Obtain the partial R and square them
partr2 <- modelEffectSizes(mlr_glob)
partr2$Effects[,4]^2



#### Figure 3: Spatial patterns in the Miellajokka catchment ####
#Read file from the spatial sampling in Miellajokka
#It can be found here: https://doi.org/10.5878/pxa2-vy55
sp <- read_csv2("spatial_miella_submission.csv") #Beware, this csv file was processed by the data service and has "," as decimal delimiter

####THIS IS TO MAKE THE K VS CO2 SPACE.
#### TO PLOT THE ISOLINES IN FIGURE 3A
pco2 <- seq(from=3800, to= 380, by=-5)
ks <- seq(from= 1, to=550, length.out = length(pco2) )
kpco2 <- list( ks, pco2)
atmCO2 <- 0.000380

#I basically do a matrix with all given values of pco2 and k and fco2
x <- matrix( nrow=length(pco2), ncol=length(pco2), dimnames = kpco2) 

#I use the mean temp of 9.2 to calculate the henrys
Henrys <- (exp(-58.0931+(90.5069*(100/(9.2+273.15)))+(22.294*log((9.2+273.15)/100))))*1000 

for(i in 1:length(pco2)){
  for( j in 1:length(ks)){
    kco2 <- (((1100/599.42)^-0.5) * ks[i]) #I use the mean temp of 9.2 to calculate the schmidt number
    x[i,j] <- kco2*Henrys*(pco2[j]/1e+6-atmCO2)*12.0107
  }
}

tx <- setNames(melt(x), c('k', 'co2', 'fco2'))
tx100<- subset(tx, fco2 < 50)

#Figure 3a
## Black and white version, export pdf 8x5 inch
ggplot() + 
  geom_contour(data=tx100, aes(y=co2, x=k,  z=fco2), binwidth = 10, linetype=3, colour="grey10")+
  geom_contour(data=tx, aes(y=co2, x=k,  z=fco2), binwidth = 50, linetype=2, colour="grey10")+
  geom_point(data = sp, aes(y=CO2_ppm, x=k600_md), size=5, alpha=.8)+
  geom_hline(aes(yintercept=380), linetype=1, color="gray30")+
  scale_y_continuous(limits = c(0, 4150), expand = c(0, 0), breaks = seq(0, 4000, by = 1000)) +
  scale_x_continuous(limits = c(0, 550),expand = c(0, 0), breaks = seq(0, 550, by = 100)) +
  labs(y=expression(paste(italic(p),CO[2], " (ppm)")), x=expression(k[600]~~(m~d^-1)))+
  theme(axis.text = element_text(size=19), axis.title = element_text(size=19), plot.subtitle = element_text(size=13),
        legend.key.height = unit(2.35, "cm"),
        legend.key.width = unit(1, "cm"))


##Color gradient version, export as pdf 8x5 inch
ggplot() + 
   geom_tile(data=tx, aes(x=k, y=co2, fill=fco2), show.legend= T) +
    scale_fill_gradientn(colours = c("white","yellow2","orange", "darkorange2", "red","firebrick4"), 
                       values = rescale(c(0,20,70,130,200,700)), guide = "colorbar", name=expression(paste( CO[2]," Flux (gC ",m^-2,d^-1,")")),
                      labels = c("0", "50", "100","200", "400", "600","800"),
                     breaks = c(0, 50, 100,200, 400,600, 800))+
  geom_contour(data=tx100, aes(y=co2, x=k,  z=fco2), binwidth = 10, linetype=3, colour="grey10")+
  geom_contour(data=tx, aes(y=co2, x=k,  z=fco2), binwidth = 50, linetype=2, colour="grey10")+
  geom_point(data = sp, aes(y=CO2_ppm, x=k600_md), size=5, alpha=.8)+
  scale_y_continuous(limits = c(0, 4150), expand = c(0, 0), breaks = seq(0, 4000, by = 1000)) +
  scale_x_continuous(limits = c(0, 550),expand = c(0, 0), breaks = seq(0, 550, by = 100)) +
  geom_hline(aes(yintercept=380), linetype=1, color="gray30")+
  labs(y=expression(paste(italic(p),CO[2], " (ppm)")), x=expression(k[600]~~(m~d^-1)))+
  theme(axis.text = element_text(size=19), axis.title = element_text(size=19), plot.subtitle = element_text(size=13),
        legend.key.height = unit(2.35, "cm"),
        legend.key.width = unit(1, "cm"))



#Figure 3b:
# Firs we subset by k threshold
wet_klow <- subset(sp, subset = k600_md <= 42)
wet_khigh <- subset(sp, subset = k600_md >42)

#Figure 3b,  export as 8x5 inch and pdf
ggplot()+
  geom_point(data=wet_klow, aes(x=wet_areas_percentage, y=CO2_ppm), size=5, alpha=0.8,  color="forestgreen")+ 
  geom_point(data=wet_khigh, aes(x=wet_areas_percentage, y=CO2_ppm), size=6, alpha=0.8,shape=18, color="darkgoldenrod3")+
  geom_smooth(data=wet_klow,aes(x=wet_areas_percentage, y=CO2_ppm),method="lm", color="forestgreen", alpha=0.2, fill="darkolivegreen3")+
  geom_smooth(data=wet_khigh,aes(x=wet_areas_percentage, y=CO2_ppm),method="lm", color ="darkgoldenrod3", alpha=0.2, fill="darkgoldenrod1")+
  scale_y_continuous(limits = c(0,4150), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 105),expand = c(0, 0), breaks = seq(0, 100, by = 25)) +
  labs(x="% wet area", y=expression("pCO"[2]*" (ppm)"))+
  theme(axis.text = element_text(size=19), axis.title = element_text(size=19), plot.subtitle = element_text(size=13),
        legend.key.height = unit(2.35, "cm"),
        legend.key.width = unit(1, "cm"))

#I add the equations and labels in another illustrator software, the lm are:
lmlow <- lm(CO2_ppm~wet_areas_percentage, data=wet_klow)
summary(lmlow)

lmhigh<- lm(CO2_ppm~wet_areas_percentage, data=wet_khigh)
summary(lmhigh)


#Figure 3c,   export as 8x5 inch and pdf
ggplot(data= sp, aes(x=k600_md, y=FluxCO2_gm2_d1, color=CO2_ppm))+
  geom_point(size=8,alpha=0.8, shape=20)+
  scale_color_gradient(low="dodgerblue", high="firebrick2")+
  scale_y_continuous(limits = c(0,62), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 550),expand = c(0, 0), breaks = seq(0, 500, by = 100)) +
  labs(colour=expression(paste(italic(p),CO[2], " (ppm)")), y=expression(CO[2]~evasion~(g~C~m^-2~d^-1)), x=expression(k[600]~~(m~d^-1)) )+
  theme(legend.position= c(0.87,0.72),legend.key.height = unit(0.9, "cm"), 
        axis.text = element_text(size=19), axis.title = element_text(size=19),
        legend.text = element_text(size=19), legend.title = element_text(size=19),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.8, "cm"))   

#Figure 3d,   export as 8x5 inch and pdf
ggplot(data= sp, aes(x=CO2_ppm, y=FluxCO2_gm2_d1))+
  geom_point(size=5, alpha=.7)+
  geom_smooth( color="black", method=lm)+
  gglm+
  scale_y_continuous(limits = c(0,62), expand = c(0, 0))+
  scale_x_continuous(limits = c(0,3600), expand = c(0, 0))+
  labs(x=expression(paste(italic(p),CO[2], " (ppm)")), y=expression(CO[2]*" evasion (g C "*m^-2*d^-1*")" ) )+
  theme(axis.text = element_text(size=19), axis.title = element_text(size=19))


#### Figure 4: Bootstraping the co2 evasion distribution ####

####This is to do the bootstrapping of the spatial distribution

x1 <-sp$FluxCO2_gm2_d1


size=168 #Maximum size
n=200 #Number of times each size
boot1 <- data.frame( n= rep(1:size, each=200), fco2=double(size*n))

#We do the sampling with replacement in this for
set.seed(2)
for(i in 1:(size*n)){
  boot1$fco2[i]<-  median( sample(x1, boot1$n[i], replace = T) )
}

#Figure 4a, the bootstrapping: export 9x6inch
ggplot()+
  geom_point(data=boot1,aes(x=n, y=fco2), alpha=0.7, size=5, color="gray1")+
  geom_point(aes(x=1, y=54.7), size=5, shape=17, color="red2")+ #evasion calculating with pco2 monitoring*average k
  geom_segment(aes(x=1,xend=size,y=median(x1),yend=median(x1)), colour="red", linetype=2, size=2)+
  scale_x_continuous(expand = c(0.02, 0),breaks= c(0,10,25,50,100,150))+
  scale_y_continuous(expand = c(0.02, 0),breaks= c(0,10,20,30,40,50,60))+
  labs(x="number of samples", y=expression("median "*CO[2]*" evasion (g C "*m^-2*d^-1*")" ))+
  theme(axis.text = element_text(size=19), axis.title = element_text(size=19))

# Figure 4b, Histogram inset,  export as 7x5 inch and pdf
ggplot() +
  geom_histogram(data = as.data.frame(x1), aes(x=x1, y=..density..), fill="gray15", color="gray25", bins = 60) +
  scale_x_continuous(expand = c(0, 0), limits=c(0,62))+
  scale_y_continuous(expand = c(0, 0))+
  geom_segment(aes(x=median(x1),xend=median(x1),y=0,yend=0.155), colour="red", size=2,linetype=2)+
  labs(y="density", x=expression(CO[2]*" evasion (g C "*m^-2*d^-1*")" ))+
  theme(axis.text = element_text(size=23), axis.title = element_text(size=23))



#### For the statistics in the log space for Miellajokka ####
splog <- sp %>% mutate(pco2=if_else(CO2_ppm<=380,380.5, CO2_ppm)) %>% #There are 2 sites with pco2<380, so I put them to 380
  transmute(siteID=siteID, fco2= log(abs(FluxCO2_gm2_d1), base=10), 
                          k600=log(kco2_md, base=10), pco2=log(pco2-380, base=10),
                          wet=log(wet_areas_percentage, base=10) ) %>% 
  filter(fco2>-500) %>%  na.omit()  %>%  #after log transform there are some -Inf values
  mutate(fco2s=scale(fco2), pco2s=scale(pco2), k600s=scale(k600))

#The multiple linear regression
mlr <- lm(fco2s~pco2s+k600s, data=splog)

#the coefficients of the MLR
summary(mlr)

#and the partial R2 for each IV
partr2 <- modelEffectSizes(mlr)
partr2$Effects[,4]^2
