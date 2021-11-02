################################################################################
################################################################################
# Silver Shiner Occupancy in 16 Mile Creek, Ontario, Canada
# R Code written by Karl A. Lamothe - karl.lamothe@dfo-mpo.gc.ca
# Great Lakes Laboratory for Fisheries and Aquatic Sciences
# The data were collected under Fisheries and Oceans Canada permitting 
# Edited 2021-11-02  R Version: 4.1.1
################################################################################
################################################################################
# load libraries
library(pacman)
p_load(unmarked)     # occupancy models
p_load(openxlsx)     # read in xlsx files
p_load(ggplot2)      # plotting
p_load(ggrepel)      # For plotting non-overlapping text
p_load(AICcmodavg)   # AIC comparisons
p_load(vegan)        # multivariate analyses
p_load(ape)          # Moran I  
p_load(patchwork)    # multiple ggplots per figure
p_load(ggcorrplot)   # calculate and plot pearson correlations

################################################################################
################################################################################
# personal ggplot theme
theme_me <- theme_bw() +
  theme(axis.title = element_text(size= 11, family= "sans", colour= "black"),
        axis.text.x = element_text(size= 10, family= "sans", colour= "black"),
        axis.text.y = element_text(size= 10, family= "sans", colour= "black"),
        legend.title = element_text(size= 10, family= "sans", colour= "black"),
        legend.text = element_text(size= 8, family= "sans", colour= "black"),
        strip.text = element_text(size= 11, family= "sans", colour= "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        plot.title = element_text(size= 11, family= "sans", colour= "black"),
        axis.ticks = element_line(colour="black"))

################################################################################
################################################################################
# load community data
# includes both juvenile and adult Silver Shiner in unique columns 
# summed after all hauls
community<-read.xlsx("SilverShiner16MileCreek2.xlsx",sheet = "Haul1+2+3")
colnames(community) 

# include only adult SS 
ad.SS.comm<-community
ad.SS.comm<-ad.SS.comm[-c(31,32)]

# include only juvenile SS
juv.SS.comm<-community
juv.SS.comm<-juv.SS.comm[-c(30,32)]

# adult and juvenile SS as unique 'species' (i.e. df column) 
ad.juve.comm<-community
ad.juve.comm<-ad.juve.comm[-32]

# full comm no separation of juveniles and adults
full.comm<-community
full.comm<-full.comm[-c(30,31)]

################################################################################
# load habitat data
################################################################################
Habitat<-read.xlsx("SilverShiner16MileCreek2.xlsx", sheet = "Habitat")
colnames(Habitat) 
str(Habitat) 

#convert year to factor
Habitat$Year<-as.character(Habitat$Year)
Habitat$Year<-as.factor(Habitat$Year)
Habitat$Type<-as.factor(Habitat$Type)

# generate means and standard deviations for habitat characteristics
aggregate(Habitat$Depth, by=list(Habitat$Year), mean, na.rm = T)
aggregate(Habitat$Depth, by=list(Habitat$Year), sd, na.rm = T)
aggregate(Habitat$WV, by=list(Habitat$Year), mean, na.rm = T)
aggregate(Habitat$WV, by=list(Habitat$Year), sd, na.rm = T)
aggregate(Habitat$Conductivity, by=list(Habitat$Year), mean, na.rm=T)
aggregate(Habitat$Conductivity, by=list(Habitat$Year), sd, na.rm=T)
aggregate(Habitat$pH, by=list(Habitat$Year), mean, na.rm=T)
aggregate(Habitat$pH, by=list(Habitat$Year), sd, na.rm=T)
aggregate(Habitat$DO, by=list(Habitat$Year), mean, na.rm=T)
aggregate(Habitat$DO, by=list(Habitat$Year), sd, na.rm=T)
aggregate(Habitat$Water_Temp, by=list(Habitat$Year), mean, na.rm=T)
aggregate(Habitat$Water_Temp, by=list(Habitat$Year), sd, na.rm=T)

aggregate(Habitat$WV, 
          by=list(Habitat$Year,Habitat$Type), mean, na.rm = T)

# Calculate correlations between variables
colnames(Habitat)
CorrVars<-cbind.data.frame(Habitat[,c(7,8,10,11,12,13)])
colnames(CorrVars)<-c("Site velocity","Depth","Water temp.",
                      "Conductivity","DO","pH")
CorrVars2<-CorrVars
CorrVars2<-cbind.data.frame(CorrVars2, Year=Habitat$Year, Type=Habitat$Type)
CorrVars<-CorrVars[complete.cases(CorrVars), ]

# Generate correlation plot with pearson correlations
correlationplot<-ggcorrplot(cor(CorrVars,method="pearson"), hc.order = TRUE, outline.col = "black",
                            type = "lower", lab=TRUE,
                            ggtheme = ggplot2::theme_gray,
                            insig = "blank")+
  theme(panel.border = element_rect(colour = "black", fill=NA),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

#plot and export
#png("Figure S1.png", width = 5, height = 5, units = 'in', res = 500)
correlationplot
#dev.off()

# Perform kruskal-wallis test on habitat characteristics
#kruskal.test(WV ~ Year, data = Habitat)
#kruskal.test(Depth ~ Year, data = Habitat)
#kruskal.test(pH ~ Year, data = Habitat, na.action = "na.omit")
#kruskal.test(DO ~ Year, data = Habitat, na.action = "na.omit")
#kruskal.test(Conductivity ~ Year, data = Habitat, na.action = "na.omit")
#kruskal.test(Water_Temp ~ Year, data = Habitat, na.action = "na.omit")

# Create boxplots for habitat variables
# Depth
Depthbox<-ggplot(Habitat, aes(y=Depth, x=as.character(Year)))+
  geom_boxplot(color="black")+
  geom_jitter(color="black", width=0.1, size=0.5) +
  labs(y="Depth (m)", x="")+
  theme_me+
  ylim(0,1.25)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Depthbox

# Water temperature
Temp<-ggplot(Habitat, aes(y=Water_Temp, x=as.character(Year)))+
  geom_boxplot(color="black")+
  geom_jitter(color="black", width=0.1, size=0.5) +
  labs(y="Water temp. (°C)", x="")+
  theme_me+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Temp

#Water velocity
WaterVel<-ggplot(Habitat, aes(y=WV, x=as.character(Year)))+
  geom_boxplot(color="black", size=0.5)+
  geom_jitter(color="black", width=0.1, size=0.5) +
  labs(y="Velocity (m/s)", x="")+
  theme_me+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
WaterVel

# Dissolved oxygen
DissO<-ggplot(Habitat, aes(y=DO, x=as.character(Year)))+
  geom_boxplot(color="black", outlier.size=0.75)+
  geom_jitter(color="black", width=0.1, size=0.5) +
  labs(y="Dissolved oxygen (mg/L)", x="Year")+
  theme_me+
  theme(legend.position = "none")
DissO

#pH
pHfigu<-ggplot(Habitat, aes(y=pH, x=as.character(Year)))+
  geom_boxplot(color="black", outlier.size=0.75)+
  geom_jitter(color="black", width=0.1, size=0.5) +
  labs(y="pH", x="")+
  theme_me+
  theme(legend.position = "none")
pHfigu

#Conductivity
Cond<-ggplot(Habitat, aes(y=Conductivity, x=as.character(Year)))+
  geom_boxplot(color="black", outlier.size=0.75)+
  geom_jitter(color="black", width=0.1, size=0.5) +
  labs(y="Conductivity (µs/cm)", x="")+
  theme_me+
  theme(legend.position = "none")
Cond

#Depth by habitat type
Dept.type.hbox<-ggplot(Habitat, aes(y=Depth, x=as.character(Type)))+
  geom_boxplot(color="black")+
  geom_jitter(color="black", width=0.1, size=0.5) +
  labs(y="Depth (m)", x="")+
  geom_point(data = data.frame(x = c(1,1,2,2,3,3), 
                               y = c(0.83,0.625,
                                     0.3075,0.2683,
                                     0.7675,0.568)),
             aes(x=x, y=y),
             size=2,
             shape = c(17,15,17,15,17,15))+
  theme_me+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Dept.type.hbox

#Water velocity by habitat type
Water.type.Vel<-ggplot(Habitat, aes(y=WV, x=as.character(Type)))+
  geom_boxplot(color="black", outlier.size=0.75)+
  geom_jitter(color="black", width=0.1, size=0.5) +
  geom_point(data = data.frame(x = c(1,1,2,2,3,3), 
                               y = c(0.12428571,0.036666,
                                     0.7025,0.4383333,
                                     0.57875,0.211)),
             aes(x=x, y=y),
             size=2,
             shape = c(17,15,17,15,17,15))+
  labs(y="Velocity (m/s)", x="Habitat type")+
  theme_me
Water.type.Vel

# Combine figures and export
#png("Figure 4.png", width = 7.5, height = 4, units = 'in', res = 800)
Depthbox + WaterVel + Temp+ Dept.type.hbox +  Cond + DissO + pHfigu +Water.type.Vel +
  plot_layout(ncol=4)
#dev.off()

################################################################################
################################################################################
# Load juvenile, adult, and combined silver shiner occupancy data format
################################################################################
################################################################################
SSdata.Adult<-read.xlsx("SilverShiner16MileCreek2.xlsx",
                        sheet = "OccurrenceAdults")
head(SSdata.Adult)

SSdata.YOY<-read.xlsx("SilverShiner16MileCreek2.xlsx",
                      sheet = "OccurrenceJuveniles")
head(SSdata.YOY)

SSdata.Comb<-read.xlsx("SilverShiner16MileCreek2.xlsx",
                      sheet = "OccurrenceComb")
head(SSdata.Comb)

################################################################################
# Testing for spatial autocorrelation in habitat variables using Moran I
################################################################################
colnames(Habitat)
Habitat2<-Habitat[c(4,5,7,8)]
Habitat2<-Habitat2[complete.cases(Habitat2), ]

# Calculate geographic distances and inverse distance matrix
SS_dists <- as.matrix(dist(cbind(Habitat2$Longitude, Habitat2$Latitude)))
SS.dists.inv <- 1/SS_dists
diag(SS.dists.inv) <- 0
SS.dists.inv[1:5, 1:5]

# Test for spatial autocorrelation
Moran.I(log(Habitat2$WV+1), SS.dists.inv, alternative="greater")
Moran.I(log(Habitat2$Depth+1), SS.dists.inv ,alternative="greater")


commtype<-cbind.data.frame(community, type=Habitat$Type)
aggregate(ad.juve.comm$Luxilus.cornutus, by=list(commtype$type),sum)
aggregate(ad.juve.comm$Luxilus.chrysocephalus, by=list(commtype$type),sum)
aggregate(community$Notropis.photogenis, by=list(commtype$type),sum)

aggregate(ad.juve.comm$SS.Adult, by=list(commtype$type),sum)
aggregate(ad.juve.comm$SS.Juvenile[ad.juve.comm$Year=="2016"], 
          by=list(commtype$type[ad.juve.comm$Year=="2016"]),sum)

98/(98+704+1202)
################################################################################
########################### Occupancy models ###################################
############################## Juveniles #######################################
################################################################################
# Calculate Naive estimate of occupancy total
Naive.SS<-sum(SSdata.YOY[3:5], na.rm = T)
Naive.SS/length(SSdata.YOY$Field.Number) #0.5217391

# Split into two years to look at sampling patterns
JVSS2011<-SSdata.YOY[SSdata.YOY$Year=="2011",]
JVSS2016<-SSdata.YOY[SSdata.YOY$Year=="2016",]

#Calculate Naive estimates of occupancy separately
Naive.SS2011<-sum(JVSS2011[3:5], na.rm = T)
Naive.SS2011/length(JVSS2011$Field.Number) #0.1666667 - rarely captured

#2016
Naive.SS2016<-sum(JVSS2016[3:5], na.rm = T)
Naive.SS2016/length(JVSS2016$Field.Number) #0.9090909 - frequently captured

################################################################################
# Collect the data for the occupancy modeling
# SS frame is the data frame of juvenile observations at each site with 3 hauls
# Habitat is the site specific covariates
################################################################################
SS.frame<-cbind.data.frame(Haul.1 = SSdata.YOY$Haul.1, 
                           Haul.2 = SSdata.YOY$Haul.2, 
                           Haul.3 = SSdata.YOY$Haul.3)
SS.frame

# Occupancy model data frame using 'unmarked'
SS.umf<-unmarkedFrameOccu(y=SS.frame, siteCovs = Habitat)

# must remove row 15 due to missing depth and water velocity measurements
SS.umf<-SS.umf[-c(15),]
head(SS.umf)
str(SS.umf)

################################################################################
# Plot observations
################################################################################
#pdf('SS.juv.detections.pdf',  width=2, height=3)
plot(SS.umf, panels=1)
#dev.off()

################################################################################
# Model building
################################################################################
fm<-list()
fm[[1]] <- occu(~1 ~1, SS.umf) #intercept model
#fm[[2]] <- occu(~WV ~1, SS.umf)  #inadequate estimates [see Standard errors]
#fm[[3]] <- occu(~Depth ~1, SS.umf) #inadequate estimates [see Standard errors]
#fm[[4]] <- occu(~Year ~1, SS.umf) #did not converge
#fm[[5]] <- occu(~Type ~1, SS.umf) #did not converge

# No Model Selection for detection
# occupancy
fm[[2]] <- occu(~1 ~Depth, SS.umf)
fm[[3]] <- occu(~1 ~WV, SS.umf)
fm[[4]] <- occu(~1 ~Year, SS.umf)

# Model Selection
Mods <- paste("fm", 1:4, sep = "")
Z<-as.data.frame(aictab(cand.set=fm, modnames=Mods, sort=T, second.ord=T, LL=T))
Z 

# report coefficients
coef(fm[[4]])
coef(fm[[3]])
coef(fm[[1]])
coef(fm[[2]])

# Chi-square test
chisq <- function(mod) {
  observed <- getY(mod@data)
  expected <- fitted(mod)
  chisq <- sum((observed - expected)^2 / expected, na.rm=T)
}

# set.seed
set.seed(4256)
fm4 <- occu(~1 ~Year, SS.umf) 
(pb <- parboot(fm4, chisq, nsim=999, report=1)) #bootstrapping

################################################################################
# Predictions
################################################################################
Predictpsi<-predict(fm4, type="state", Habitat)
Predictp<-predict(fm4, type="det", Habitat)
head(Predictp)

# Averages
mean(Predictp$Predicted) #0.8852966 
sd(Predictp$Predicted)/sqrt(length(Predictp$Predicted)) #0
mean(Predictpsi$Predicted) #0.5263106
sd(Predictpsi$Predicted)/sqrt(length(Predictpsi$Predicted)) # 0.05482644

# probability of detecting juvenile silver shiners at least once after 3 surveys
# given that the site is occupied
1-(1-0.8852966)^3

#################################################################################
# Covariate prediction
################################################################################
Year<-c(2011,2016)
newdata<-cbind.data.frame(Year=as.factor(Year))
occu<-predict(fm4, type='state', newdata=newdata)
occu<-cbind.data.frame(occu,Year=Year)
str(occu)

# convert to factor
occu$Year<-as.character(occu$Year)
occu$Year<-as.factor(occu$Year)

# plot occupancy as a function of year
occ.plt1<-ggplot(occu, aes(x=Year, y=Predicted)) + 
  ylim(0,1)+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0) +
  geom_point() + ylab("Occupancy probability") + xlab("Year") + theme_me +
  ggtitle("Juveniles")
occ.plt1

################################################################################
########################### Occupancy models ###################################
################################ Adults ########################################
################################################################################
# Calculate Naive estimate of occupancy
Naive.SS<-sum(SSdata.Adult[3:5], na.rm = T)
Naive.SS/length(SSdata.Adult$Field.Number) #0.4347826

#Split into two years
ADSS2011<-SSdata.Adult[SSdata.Adult$Year=="2011",]
ADSS2016<-SSdata.Adult[SSdata.Adult$Year=="2016",]

#Calculate Naive estimates of occupancy seperately
Naive.SS2011<-sum(ADSS2011[3:5], na.rm = T)
Naive.SS2011/length(ADSS2011$Field.Number) #0.5

#2016
Naive.SS2016<-sum(ADSS2016[3:5], na.rm = T)
Naive.SS2016/length(ADSS2016$Field.Number) #0.3636364

################################################################################
# Making the response variable data frame for the adult occupancy model
################################################################################
SS.frame2<-cbind.data.frame(Haul.1 = SSdata.Adult$Haul.1, 
                            Haul.2 = SSdata.Adult$Haul.2, 
                            Haul.3 = SSdata.Adult$Haul.3)

# Occupancy model data frame using 'unmarked'
SS.umf2<-unmarkedFrameOccu(y=SS.frame2, siteCovs = Habitat)
SS.umf2<-SS.umf2[-c(15),] # as done above
head(SS.umf2)
str(SS.umf2)

################################################################################
# Plot observations
################################################################################
#pdf('SS.ad.detections.pdf',  width=2, height=3)
plot(SS.umf2, panels=1)
#dev.off()

################################################################################
# Model building
################################################################################
fm<-list()
fm[[1]] <- occu(~1 ~1, SS.umf2) #intercept model
#fm[[2]] <- occu(~WV ~1, SS.umf2) # inadequate estimates [Standard errors]
#fm[[3]] <- occu(~Depth ~1, SS.umf2)# inadequate estimates [Standard errors]
fm[[2]] <- occu(~Year ~1, SS.umf2)

# Model Selection
Mods <- paste("fm", 1:2, sep = "")
Z<-as.data.frame(aictab(cand.set=fm, modnames=Mods, sort=T, second.ord=T, LL=T))
Z

#
fm[[3]] <- occu(~1 ~Depth, SS.umf2)
fm[[4]] <- occu(~1 ~WV, SS.umf2)
fm[[5]] <- occu(~1 ~Year, SS.umf2)

# Model Selection
Mods <- paste("fm", 1:5, sep = "")
Z<-as.data.frame(aictab(cand.set=fm, modnames=Mods, sort=T, second.ord=T, LL=T))
Z

summary(fm[[3]]) 

coef(fm[[3]])
coef(fm[[4]])
coef(fm[[1]])
coef(fm[[5]])
coef(fm[[2]])

fm3 <- occu(~ 1 ~Depth, SS.umf2) 
(pb <- parboot(fm3, chisq, nsim=999, report=1))

################################################################################
#predictions
Predictpsi<-predict(fm3, type="state", Habitat2)
Predictp<-predict(fm3, type="det", Habitat2)
head(Predictp)

#Averages
mean(Predictp$Predicted) #0.8118967 
sd(Predictp$Predicted)/sqrt(length(Predictp$Predicted)) #0
mean(Predictpsi$Predicted) #0.4251981
sd(Predictpsi$Predicted)/sqrt(length(Predictpsi$Predicted)) #0.03345092

# probability of detecting adult silver shiners at least once after 3 surveys
# given that the site is occupied
1-(1-0.8118967)^3

#################################################################################
# Covariate plot
Depth<-seq(from=min(Habitat2$Depth,na.rm=TRUE),
        to=max(Habitat2$Depth,na.rm=TRUE), by=0.01)
newdata<-as.data.frame(Depth)
occu<-predict(fm3, type='state', newdata=newdata)

occ.plt2<-ggplot(occu)+
  geom_ribbon(aes(x=Depth,ymin=lower, ymax=upper), fill="grey75")+
  geom_line(aes(x=Depth, y=Predicted), size=0.5)+
  xlab("Depth (m)")+
  ylim(0,1)+ 
  ylab("")+
  ggtitle("Adults")+
  scale_x_continuous(limits=c(0,1.25),
                     breaks=c(0.00,0.25,0.50,0.75,1.00,1.25),# next line ensures 2 decimals on x axis
                     labels=function(x){sprintf("%.2f", x)})+
  theme_me+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
occ.plt2

################################################################################
########################### Occupancy models ###################################
############################### Combined #######################################
################################################################################
# Calculate Naive estimate of occupancy
Naive.SS<-sum(SSdata.Comb[3:5], na.rm = T)
Naive.SS/length(SSdata.Comb$Field.Number) #0.6956522

#Split into two years
ADSS2011<-SSdata.Comb[SSdata.Comb$Year=="2011",]
ADSS2016<-SSdata.Comb[SSdata.Comb$Year=="2016",]

#Calculate Naive estimates of occupancy separately
Naive.SS2011<-sum(ADSS2011[3:5], na.rm = T)
Naive.SS2011/length(ADSS2011$Field.Number) #0.5

#2016
Naive.SS2016<-sum(ADSS2016[3:5], na.rm = T)
Naive.SS2016/length(ADSS2016$Field.Number) #0.9090909

################################################################################
# Making the response variable data frame for the combined occupancy model
################################################################################
SS.frame3<-cbind.data.frame(Haul.1 = SSdata.Comb$Haul.1, 
                            Haul.2 = SSdata.Comb$Haul.2, 
                            Haul.3 = SSdata.Comb$Haul.3)

# Occupancy model data frame using 'unmarked'
SS.umf3<-unmarkedFrameOccu(y=SS.frame3, siteCovs = Habitat)
head(SS.umf3)
SS.umf3<-SS.umf3[-c(15),]
str(SS.umf3)

################################################################################
# Plot observations
################################################################################
#pdf('SS.comb.detections.pdf',  width=2, height=3)
plot(SS.umf3, panels=1)
#dev.off()

################################################################################
# Model building
################################################################################
fm<-list()
fm[[1]] <- occu(~1 ~1, SS.umf3) #intercept model
#fm[[2]] <- occu(~WV ~1, SS.umf3) # inadequate estimates [Standard errors]
fm[[2]] <- occu(~Depth ~1, SS.umf3)# 
fm[[3]] <- occu(~Year ~1, SS.umf3)# 

Mods <- paste("fm", 1:3, sep = "")
Z<-as.data.frame(aictab(cand.set=fm, modnames=Mods, sort=T, 
                        second.ord=T, LL=T))
Z

# 
fm[[4]] <- occu(~1 ~Depth, SS.umf3)
fm[[5]] <- occu(~1 ~WV, SS.umf3)
fm[[6]] <- occu(~1 ~Year, SS.umf2)

# Model Selection
Mods <- paste("fm", 1:6, sep = "")
Z<-as.data.frame(aictab(cand.set=fm, modnames=Mods, sort=T, second.ord=T, LL=T))
Z

summary(fm[[5]]) 

coef(fm[[5]])
coef(fm[[6]])
coef(fm[[1]])
coef(fm[[4]])
coef(fm[[3]])
coef(fm[[2]])

fm5 <- occu(~ 1 ~WV, SS.umf3) 
(pb <- parboot(fm5, chisq, nsim=999, report=1))

################################################################################
#predictions
Predictpsi<-predict(fm5, type="state", Habitat2)
Predictp<-predict(fm5, type="det", Habitat2)
head(Predictp)

#Averages
mean(Predictp$Predicted) #0.8534134
sd(Predictp$Predicted)/sqrt(length(Predictp$Predicted)) #0
mean(Predictpsi$Predicted) #0.691256
sd(Predictpsi$Predicted)/sqrt(length(Predictpsi$Predicted)) #0.02715894

# probability of detecting adult silver shiners at least once after 3 surveys
# given that the site is occupied
1-(1-0.8534134)^3

#################################################################################
# Covariate plot
WV<-seq(from=min(Habitat$WV,na.rm=TRUE),
        to=max(Habitat$WV,na.rm=TRUE), by=0.01)
newdata<-as.data.frame(WV)
occu<-predict(fm5, type='state', newdata=newdata)

occ.plt3<-ggplot(occu)+
  geom_ribbon(aes(x=WV,ymin=lower, ymax=upper), fill="grey75")+
  geom_line(aes(x=WV, y=Predicted), size=0.5)+
  xlab("Water velocity (m/s)")+
  ylim(0,1)+ ylab("")+
  ggtitle("Combined")+
  theme_me+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
occ.plt3

#################################################################################
# plot the occupancy covariate plots together
#################################################################################
Figure3<-occ.plt1+occ.plt2 + occ.plt3+ 
  plot_layout(widths = c(1, 2,2))
 
#png("Figure 3.png", width = 6, height = 2.5, units = 'in', res = 800)
Figure3 # patchwork makes this really easy to do!
#dev.off()

################################################################################
################################################################################
############################## Community Analysis ##############################
################################################################################
################################################################################
# Make Data frames
Fish.counts2<-as.data.frame(t(aggregate(full.comm[3:30],
                                     by=list(Habitat$Year),sum)))
head(Fish.counts2)
Fish.counts2<-as.data.frame(Fish.counts2[-1,]) #remove top row
colnames(Fish.counts2)<-c("2011","2016") #revise column names

# counts for ggplot
Fish.countsgg<-cbind.data.frame(Species = c(rownames(Fish.counts2),
                                            rownames(Fish.counts2)),
                                   Count = c(Fish.counts2$`2011`, 
                                             Fish.counts2$`2016`),
                                   Year = c(rep("2011",28),
                                            rep("2016",28)))
head(Fish.countsgg)
str(Fish.countsgg)
Fish.countsgg$Year<-as.factor(Fish.countsgg$Year)
Fish.countsgg$Count<-as.numeric(Fish.countsgg$Count)

# replace period with a space in the species names
Fish.countsgg$Species<-gsub(".", " ", Fish.countsgg$Species, fixed=TRUE)

# Total number of fish caught
sum(Fish.countsgg$Count[Fish.countsgg$Year=="2011"])
sum(Fish.countsgg$Count[Fish.countsgg$Year=="2016"])

# Total number of individuals per species caught
aggregate(Fish.countsgg$Count, by=list(Fish.countsgg$Species), sum)

# Convert this to relative abundance per year
Counts2011<-Fish.countsgg$Count[Fish.countsgg$Year=="2011"]
Counts2016<-Fish.countsgg$Count[Fish.countsgg$Year=="2016"]
Relative2011<-Counts2011/sum(Counts2011)
Relative2016<-Counts2016/sum(Counts2016)
Relative<-c(Relative2011,Relative2016)

# Combine the data frames
Fish.countsgg<-cbind(Fish.countsgg, Relative)

# Only take the most prevalent species
Fish.countsgg.red<-Fish.countsgg[c(1,3,4,6,8,12,13,14,17,19,23,26,27,28,
                                   29,31,32,34,36,40,41,42,45,47,51,54,55,56),]

# Plot
FigureS2<-ggplot(Fish.countsgg.red, aes(reorder(Species, Relative), 
                              fill=Year, Relative))+ 
  geom_bar(stat="identity", position="dodge")+
  ylab("Proportion of annual catch") + 
  xlab("Species") +
  coord_flip()+
  scale_fill_manual(values=c("black","gray"))+
  scale_y_continuous(limits=c(0,0.5), labels=function(x){sprintf("%.2f", x)})+
  theme_me+
  theme(axis.text.y = element_text(face="italic"),
        legend.position = c(.9, .3),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_blank())

# Export figure
#png("Figure S2.png", width = 6, height = 3, units = 'in', res = 800)
FigureS2
#dev.off()

# create capture summary of silver shiner, common shiner, and striped shiner
aggregate(ad.juve.comm$SS.Adult,by=list(ad.juve.comm$Year, Habitat$Type),sum)
aggregate(ad.juve.comm$Luxilus.cornutus,by=list(ad.juve.comm$Year, Habitat$Type),sum)
aggregate(ad.juve.comm$Luxilus.chrysocephalus,by=list(ad.juve.comm$Year, Habitat$Type),sum)
aggregate(ad.juve.comm$Notropis.rubellus,by=list(ad.juve.comm$Year, Habitat$Type),sum)

# Silver shiner captures
SSpres<-data.frame(Age=c(rep("Juvenile SS",6),rep("Adult SS",6)),
                   Year=c("2011","2016","2011","2016","2011","2016",
                          "2011","2016","2011","2016","2011","2016"),
                   Habitat=c("Pool","Pool","Riffle","Riffle","Run","Run",
                             "Pool","Pool","Riffle","Riffle","Run","Run"),
                   Count=c(4,380,0,93,2,898,200,120,3,2,224,78))

# Common shiner, striped shiners, rosyface shiners
Cooccur<-data.frame(Species=c(
                              rep("Common Shiner",6),rep("Striped Shiner",6),
                              rep("Rosyface Shiner",6)),
                    Year=c(
                           "2011","2016","2011","2016","2011","2016",
                           "2011","2016","2011","2016","2011","2016",
                           "2011","2016","2011","2016","2011","2016"),
                    Habitat=c(
                              "Pool","Pool","Riffle","Riffle","Run","Run",
                              "Pool","Pool","Riffle","Riffle","Run","Run",
                              "Pool","Pool","Riffle","Riffle","Run","Run"),
                    Count=c(
                            348,1250,8,404,270,1551,173,17,14,11,196,227,
                            34,6,0,99,9,163))

# silver shiner capture
count.type<-ggplot(SSpres, aes(x=Habitat,y=Count, fill=Age))+
  geom_bar(stat="identity")+
  facet_wrap(~Year, scales = "free")+
  scale_fill_manual(name="Life-stage",
                    values=c("darkgrey","lightgrey"))+
  theme_me+
  xlab("")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
count.type

# cooccurring species capture
count.type2<-ggplot(Cooccur, aes(x=Habitat,y=Count, fill=Species))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values=c("#0072B2","#CC79A7",
                             "#F0E442"))+
  facet_wrap(~Year, scales = "free_y")+
  theme_me+
  theme(strip.text.x = element_blank())
count.type2

# Export figure
#png("Figure XX.png", width = 6, height = 4, units = 'in', res = 800)
count.type/ count.type2+
  plot_layout(guides = "collect")
#dev.off()

################################################################################
################################################################################
## Redundancy Analysis
################################################################################
################################################################################
# transform data
head(ad.SS.comm)

# hellinger trasformation
ad.SS.comm.trans<-decostand(ad.SS.comm[3:30], method = "hellinger")
juv.SS.comm.trans<-decostand(juv.SS.comm[3:30], method = "hellinger")
ad.juve.comm.trans<-decostand(ad.juve.comm[3:31], method = "hellinger")
total.comm.trans<-decostand(full.comm[3:30], method = "hellinger")

#choose variables for RDA
head(Habitat)
Habitat.RDA<-cbind.data.frame(Habitat[c(6,7,8)])
str(Habitat.RDA)

# remove row with missing data
Habitat.RDA<-Habitat.RDA[-c(15),]
ad.juve.comm.trans<-ad.juve.comm.trans[-c(15),]
################################################################################
# Run the model for adults and juveniles separated
################################################################################
rda.fin <-rda(ad.juve.comm.trans ~ WV + Depth, data = Habitat.RDA)
summary(rda.fin)

vif.cca(rda.fin) #variable inflation factor

# Calculate adjusted R2 of RDA
RsquareAdj(rda.fin)$adj.r.squared

# Proportion of variance explained per component
rda.fin$CCA$eig/sum(rda.fin$CCA$eig)

# screeplot
screeplot(rda.fin)

###################################### Triplot #################################
# site scores
scores <- data.frame(Habitat.RDA,rda.fin$CCA$u)

# Species scores
vscores <- data.frame(rda.fin$CCA$v)

# Covariate scores
var.score <- data.frame(rda.fin$CCA$biplot) 
var.score$variables <- c("Velocity","Depth")

# Only plotting the most abundant species for visualization purposes
rownames(vscores)
colSums(ad.juve.comm[3:31])
vscores2 <- vscores[c(1,3,4,6,8,12:14,17,19,23,26:29),]
rownames(vscores2)<-c("RoBa","WhSu","BrSt","RaDa", "JoDa", "StSh","CoSh",
                      "SMB","RiCh","RoSh","FaMi","LoDa","CrCh","SSa","SSj")
cols<-c("2011" = "black", 
        "2016" = "grey") #create colour scheme for plot

# Create plot of model
ggRDA <- ggplot(scores, aes(x = RDA1, y = RDA2)) +
  geom_hline(yintercept=0,linetype="dashed",col="black")+
  geom_vline(xintercept=0,linetype="dashed",col="black")+
  geom_point(aes(col=Habitat.RDA$Year)) +
  scale_color_manual(values=cols,
                     name="Year") +
  geom_segment(data = vscores2, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow=arrow(length=unit(0.2,"cm")), 
               color = "black",inherit.aes = FALSE,lwd=0.25) +
  geom_segment(data = var.score, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length=unit(0.2,"cm")), 
               color = 'black', inherit.aes = FALSE, lwd=0.25) +
  geom_text_repel(data = vscores2, 
                  aes(x = RDA1, y = RDA2, label = rownames(vscores2)), 
                  max.overlaps = 20, size=3)+ 
  coord_fixed()+
  geom_text_repel(data = var.score, 
                  aes(x = RDA1, y = RDA2, label = variables),
                  nudge_y = ifelse(var.score$RDA2 > 0, .05, -.05),
                  segment.alpha = 0, size=3) +
  labs(x = "RDA Axis 1 (72.2%)", y = "RDA Axis 2 (27.8%)") +
  theme_me+
  theme(legend.position = c(0.85,0.25),
        legend.background = element_blank())
  
ggRDA2

################################################################################
# Run the model for adults and juveniles combined
################################################################################
total.comm.trans<-total.comm.trans[-c(15),]
rda.fin.1 <-rda(total.comm.trans ~ Depth + WV, 
                data = Habitat.RDA)
summary(rda.fin.1)

vif.cca(rda.fin.1) #variable inflation factor looks good

# Calculate adjusted R2 of RDA
RsquareAdj(rda.fin.1)$adj.r.squared

# Proportion of variance explained per component
rda.fin.1$CCA$eig/sum(rda.fin.1$CCA$eig)

# screeplot
screeplot(rda.fin.1)

# variable coefficients
coef(rda.fin.1)

###################################### Triplot #################################
# site scores
scores.1 <- data.frame(Habitat.RDA,rda.fin.1$CCA$u)

# Species scores
vscores.1 <- data.frame(rda.fin.1$CCA$v)

# Covariate scores
var.score.1 <- data.frame(rda.fin.1$CCA$biplot) 
var.score.1$variables <- c("Depth","Velocity")

# Only plotting the most abundant species for visualization purposes
rownames(vscores.1)
colSums(full.comm[3:30])
vscores2.1 <- vscores.1[c(1,3,4,6,8,12:14,17,19,23,26:28),]
rownames(vscores2.1)<-c("RoBa","WhSu","BrSt","RaDa", "JoDa", "StSh","CoSh",
                      "SMB","RiCh","RoSh","FaMi","LoDa","CrCh","SS")

# Create plot of model
ggRDA2 <- ggplot(scores.1, aes(x = RDA1, y = RDA2)) +
  geom_hline(yintercept=0,linetype="dashed",col="black")+
  geom_vline(xintercept=0,linetype="dashed",col="black")+
  geom_point(aes(col=Habitat.RDA$Year)) +
  scale_color_manual(values=cols,
                     name="Year")+
  geom_segment(data = vscores2.1, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow=arrow(length=unit(0.2,"cm")), 
               color = "black",inherit.aes = FALSE,lwd=0.25) +
  geom_segment(data = var.score.1, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length=unit(0.2,"cm")), 
               color = 'black', inherit.aes = FALSE, lwd=0.25) +
  geom_text_repel(data = vscores2.1, 
                  aes(x = RDA1, y = RDA2, 
                      label = rownames(vscores2.1)), 
                  max.overlaps = 20, size=3)+
  geom_text_repel(data = var.score.1, 
                  aes(x = RDA1, y = RDA2, label = variables),
                  nudge_y = ifelse(var.score.1$RDA2 > 0, .05, -.05),
                  segment.alpha = 0, size=3) + 
  labs(x = "RDA Axis 1 (84.8%)", y = "RDA Axis 2 (15.2%)") + 
  coord_fixed()+
  theme_me +
  theme(legend.position = "none")

ggRDA2

# Export the figure
# Using the patchwork package to align the two ggRDAs
#png("Figure 6.png", width = 4, height = 6, units = 'in', res = 800)
ggRDA / ggRDA2+
  plot_annotation(tag_levels = "A")
#dev.off()

################################################################################
################################################################################
# Community PCoA
################################################################################
################################################################################
# read in haul specific community data
Hauls1<-read.xlsx("SilverShiner16MileCreek.xlsx",  sheet = "Haul1")
Hauls12<-read.xlsx("SilverShiner16MileCreek.xlsx", sheet = "Haul1+2")
Hauls3<-read.xlsx("SilverShiner16MileCreek.xlsx",  sheet = "Haul1+2+3")
colnames(Hauls12)

Hauls1  <- Hauls1[3:32] #remove site information
Hauls12 <- Hauls12[3:32]
Hauls3  <- Hauls3[3:32]

# remove columns with zero individuals caught
# Haul 1
i <- (colSums(Hauls1, na.rm=T) != 0)
Hauls1 <- Hauls1[, i] # keep all the non-zero columns

# Haul2
i <- (colSums(Hauls12, na.rm=T) != 0)
Hauls12 <- Hauls12[, i] 

#Haul3
i <- (colSums(Hauls3, na.rm=T) != 0)
Hauls3 <- Hauls3[, i]
#########################

# Remove juvenile and adult differentiation columns
Hauls1  <-Hauls1[-c(25:26)]
colnames(Hauls1)
Hauls12 <-Hauls12[-c(26:27)]
colnames(Hauls12)
Hauls3  <-Hauls3[-c(28:29)]
colnames(Hauls3)

# Transform data for performing PCoA (relative abundance)
Hauls1.trans  <-decostand(Hauls1, method = "hellinger")
Hauls12.trans <-decostand(Hauls12, method = "hellinger")
Hauls3.trans  <-decostand(Hauls3, method = "hellinger")

# Perform PCAs
Hauls1.PCA   <- rda(Hauls1.trans)
Hauls12.PCA  <- rda(Hauls12.trans)
Hauls3.PCA   <- rda(Hauls3.trans)

# biplots
plot(Hauls1.PCA)
plot(Hauls12.PCA)
plot(Hauls3.PCA)

# screeplots
screeplot(Hauls1.PCA)
screeplot(Hauls12.PCA)
screeplot(Hauls3.PCA)

# Haul 1 and 2 Procrustes
set.seed(4256)
Haul12.pro <- procrustes(X = Hauls12.PCA, Y = Hauls1.PCA, 
                         symmetric = TRUE, scale = FALSE)
Haul12.pro

# Haul 2 and 3 Procrustes
Haul23.pro <- procrustes(X = Hauls3.PCA, Y = Hauls12.PCA, 
                         symmetric = TRUE, scale = FALSE)
Haul23.pro

# Haul 1 and 3 Procrustes
Haul13.pro <- procrustes(X = Hauls3.PCA, Y = Hauls1.PCA, 
                         symmetric = TRUE, scale = FALSE)
Haul13.pro

# Plots
par(mar=c(5,5,2,1))

# Hauls 1 and 2
plot(Haul12.pro, kind=1, type="text")
plot(Haul12.pro, kind=1, pch=20, main="", las=1, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
plot(Haul12.pro, kind=2, main="",las=1,ylim=c(0,0.15))

# Hauls 2 and 3
plot(Haul23.pro, kind=1, pch=20, main="", las=1, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), ar.col="red")
plot(Haul23.pro, kind=1, type="text")
plot(Haul23.pro, kind=2, main="", las=1, ylim=c(0,0.15))

# Hauls 2 and 3
plot(Haul13.pro, kind=1, pch=20, main="", las=1, xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), ar.col="red")
plot(Haul13.pro, kind=1, type="text")
plot(Haul13.pro, kind=2, main="", las=1, ylim=c(0,0.15))

#Test for significant correlations 
set.seed(4256)
protest(X = Hauls12.PCA, Y = Hauls1.PCA, scores="sites", permutations=9999)
protest(X = Hauls3.PCA, Y = Hauls12.PCA, scores="sites", permutations=9999)
protest(X = Hauls3.PCA, Y = Hauls1.PCA, scores="sites", permutations=9999)

# check out residuals for final procrustes
Haul3pro <- procrustes(X = Hauls3.PCA, Y = Hauls1.PCA, 
                       symmetric = TRUE, scale = FALSE)
residuals(Haul3pro)

order(residuals(Haul3pro))
plot(Haul3pro, kind=2, main="", las=1)

# procrustes plot
Haul3pro1.df<-cbind.data.frame(Haul3pro$X)
Haul3pro2.df<-cbind.data.frame(Haul3pro$Yrot)
colnames(Haul3pro2.df)<-c( "PC1","PC2")
Haul3progg.df<-rbind.data.frame(Haul3pro1.df,Haul3pro2.df)
Haul3progg.df<-cbind(Haul3progg.df, 
                     Site=rep(seq(1:46),2), 
                     ordination = c(rep("Ord3",46),
                                    rep("Ord1",46)),
                     Outliers=c(rep("black",15),"red","red",
                                rep("black",4),"red",
                                rep("black",15),"red","red",
                                rep("black",7),
                                rep("black",15),"red","red",
                                rep("black",4),"red",
                                rep("black",15),"red","red",
                                rep("black",7)))

#ggplot
# Include sites
Haul3progg.df[Haul3progg.df$Site=="16",]
Haul3progg.df[Haul3progg.df$Site=="17",]
Haul3progg.df[Haul3progg.df$Site=="22",]
Haul3progg.df[Haul3progg.df$Site=="38",]
Haul3progg.df[Haul3progg.df$Site=="39",]

Procrustesgg<-ggplot(Haul3progg.df)+
  geom_hline(yintercept=0,linetype="dashed",col="black")+
  geom_vline(xintercept=0,linetype="dashed",col="black")+
  geom_line(aes(x=PC1, y=PC2, group = Site, col=Outliers))+
  scale_colour_manual(values=c("black","red"))+
  geom_point(aes(x=PC1, y=PC2, fill=ordination), shape=21)+
  scale_fill_manual(values=c("black","gray"))+
  theme_me+ 
  coord_fixed()+
  annotate("text",x=-0.05,y=-0.09,label="39")+
  annotate("text",x=0.17,y=-0.07,label="16")+
  annotate("text",x=0.05,y=-0.18,label="22")+
  annotate("text",x=0.11,y=-0.08,label="38")+
  annotate("text",x=-0.135,y=0.03,label="17")+
  xlab("Dimension 1")+ylab("Dimension 2")+
  theme(legend.position = "none")
Procrustesgg

# residuals plot
Residuals.proc<-as.data.frame(cbind(SS=residuals(Haul3pro),
                                    Site=seq(1:46),
                                    Sitecol=c(rep("black",15),"red","red",
                                              rep("black",4),"red",
                                              rep("black",15),"red","red",
                                              rep("black",7))))

head(Residuals.proc)
str(Residuals.proc)
Residuals.proc$SS<-as.character(Residuals.proc$SS)
Residuals.proc$SS<-as.numeric(Residuals.proc$SS)
Residuals.proc$Site<-as.character(Residuals.proc$Site)
Residuals.proc$Site<-as.numeric(Residuals.proc$Site)
quantile(residuals(Haul3pro))

# ggplot
Residplotgg<-ggplot(Residuals.proc)+
  geom_hline(yintercept=0.023728182   ,linetype="dashed",col="black")+
  geom_hline(yintercept=0.037678927   ,col="black")+
  geom_hline(yintercept=0.087691908   ,linetype="dashed",col="black")+
  geom_col(aes(y=SS, x=Site, fill=Sitecol), width = 0.5)+
  scale_fill_manual(values=c("black","red"))+
  geom_vline(xintercept=24.5, color = "black")+
  scale_y_continuous(limits=c(0,0.2), labels=function(x){sprintf("%.2f", x)})+
  ylab("Residual deviance")+
  coord_flip()+
  annotate("text",x=15, y=0.175, label="2011")+
  annotate("text",x=30, y=0.175, label="2016")+
  theme_me+
  theme(legend.position = "none")
Residplotgg

# Plot the ordination plot and residual plot together
Figure5<-Procrustesgg+Residplotgg+
  plot_annotation(tag_levels = 'A')+
  plot_layout(widths = c(2, 1))
Figure5 

# Export
#png("Figure 5.png", width = 7.5, height = 4, units = 'in', res = 800)
#tiff("Figure 5.tiff", width = 7.5, height = 4, units = 'in', res = 1000)
Figure5
#dev.off()
