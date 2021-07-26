################################################################################
################################################################################
# Silver Shiner Occupancy in 16 Mile Creek, Ontario, Canada
# R Code written by Karl A. Lamothe - karl.lamothe@dfo-mpo.gc.ca
# Great Lakes Laboratory for Fisheries and Aquatic Sciences
# The data were collected under Fisheries and Oceans Canada permitting 
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
p_load(corrplot)     # calculate and plot pearson correlations

################################################################################
# personal ggplot theme
################################################################################
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
# load community data
################################################################################
# includes both juvenile and adult Silver Shiner in unique columns 
# summed after all hauls
community<-read.xlsx("SilverShiner16MileCreek.xlsx",sheet = "Haul1+2+3")
colnames(community) 

#Notropis photogenis is combined Adult and Juvenile Silver Shiner
#SS.Adult = Adult Silver Shiner
#SS.Juvenile = Juvenile Silver Shiner

# include only adult SS 
ad.SS.comm<-community
ad.SS.comm<-ad.SS.comm[-c(33,34)]

# include only juvenile SS
juv.SS.comm<-community
juv.SS.comm<-juv.SS.comm[-c(32,34)]

# adult and juvenile SS as unique 'species' (i.e. df column) 
ad.juve.comm<-community
ad.juve.comm<-ad.juve.comm[-34]

# full comm no separation of juveniles and adults
full.comm<-community
full.comm<-full.comm[-c(32,33)]

################################################################################
# load habitat data
################################################################################
Habitat<-read.xlsx("SilverShiner16MileCreek.xlsx", sheet = "Habitat")
colnames(Habitat) # SubND = substrate not determined
head(Habitat) 

# convert percent to proportion
Habitat<-cbind.data.frame(Habitat[,c(1,2,4,5,6,17,18)],
                          Habitat[,c(7:16,19,20)]/100)
head(Habitat)

################################################################################
# Load juvenile, adult, and combined silver shiner occupancy data format
################################################################################
SSdata.Adult<-read.xlsx("SilverShiner16MileCreek.xlsx",
                        sheet = "OccurrenceAdults")
head(SSdata.Adult)

SSdata.YOY<-read.xlsx("SilverShiner16MileCreek.xlsx",
                      sheet = "OccurrenceJuveniles")
head(SSdata.YOY)

SSdata.Comb<-read.xlsx("SilverShiner16MileCreek.xlsx",
                      sheet = "OccurrenceComb")
head(SSdata.Comb)

################################################################################
# Testing for spatial autocorrelation in habitat variables using Moran I
################################################################################
# Calculate geographic distances and inverse distance matrix
SS_dists <- as.matrix(dist(cbind(Habitat$Longitude, Habitat$Latitude)))
SS.dists.inv <- 1/SS_dists
diag(SS.dists.inv) <- 0
SS.dists.inv[1:5, 1:5]

# Test for spatial autocorrelation
Moran.I(log(Habitat$WV+1), SS.dists.inv, alternative="greater")
Moran.I(log(Habitat$Depth+1), SS.dists.inv ,alternative="greater")
Moran.I(log(Habitat$Cobble+1), SS.dists.inv, alternative="greater")
Moran.I(log(Habitat$Bedrock+1), SS.dists.inv ,alternative="greater")
Moran.I(log(Habitat$Gravel+1), SS.dists.inv ,alternative="greater")
Moran.I(log(Habitat$Submerged+1), SS.dists.inv ,alternative="greater")
Moran.I(log(Habitat$OpenWater+1), SS.dists.inv ,alternative="greater")

################################################################################
########################### Occupancy models ###################################
############################## Juveniles #######################################
################################################################################
# Calculate Naive estimate of occupancy total
Naive.SS<-sum(SSdata.YOY[3:5], na.rm = T)
Naive.SS/length(SSdata.YOY$Field.Number) #0.5348837

# Split into two years to look at sampling patterns
JVSS2011<-SSdata.YOY[SSdata.YOY$Year=="2011",]
JVSS2016<-SSdata.YOY[SSdata.YOY$Year=="2016",]

#Calculate Naive estimates of occupancy separately
Naive.SS2011<-sum(JVSS2011[3:5], na.rm = T)
Naive.SS2011/length(JVSS2011$Field.Number) #0.1818182 - rarely captured

#2016
Naive.SS2016<-sum(JVSS2016[3:5], na.rm = T)
Naive.SS2016/length(JVSS2016$Field.Number) #0.9047619 - frequently captured

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
# note that the warning message is OK - its just converting the variables
# in the Habitat data frame that contain characters into factors. We will
# not be using these variables in our analysis

head(SS.umf)

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
#fm[[4]] <- occu(~Submerged ~1, SS.umf) #inadequate estimates [see Standard errors]

# Because the detection models all failed to provide any useful information,
# we move forward with our modeling using an intercept model for detection
fm[[2]] <- occu(~1 ~Depth, SS.umf)
fm[[3]] <- occu(~1 ~WV, SS.umf)
#fm[[4]] <- occu(~1 ~Submerged, SS.umf) #inadequate estimates [see Standard errors]

# Model Selection
Mods <- paste("fm", 1:3, sep = "")
Z<-as.data.frame(aictab(cand.set=fm, modnames=Mods, sort=T, second.ord=T, LL=T))
Z 

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
fm3 <- occu(~1 ~WV, SS.umf) 
(pb <- parboot(fm3, chisq, nsim=999, report=1)) #bootstrapping

################################################################################
# Predictions
################################################################################
Predictpsi<-predict(fm3, type="state", Habitat)
Predictp<-predict(fm3, type="det", Habitat)
head(Predictp)

# Averages
mean(Predictp$Predicted) #0.8803756
sd(Predictp$Predicted)/sqrt(length(Predictp$Predicted)) #0
mean(Predictpsi$Predicted) #0.5358629
sd(Predictpsi$Predicted)/sqrt(length(Predictpsi$Predicted)) #0.02788745

# probability of detecting juvenile silver shiners at least once after 3 surveys
# given that the site is occupied
1-(1-0.8803756)^3

#################################################################################
# Covariate plot
################################################################################
WV<-seq(from=min(Habitat$WV,na.rm=TRUE),
        to=max(Habitat$WV,na.rm=TRUE), by=0.01)
newdata<-as.data.frame(WV)
occu<-predict(fm3, type='state', newdata=newdata)

occ.plt<-ggplot()+
  geom_ribbon(data=occu, aes(x=WV,ymin=lower, # confidence intervals
                             ymax=upper), fill="grey75")+
  geom_line(data=occu, aes(x=WV, y=Predicted), size=0.5)+
  ylim(0,1)+
  ggtitle("Juveniles")+
  xlab("Water velocity (m/s)")+
  ylab("Occupancy Probability")+
  theme_me
occ.plt

################################################################################
########################### Occupancy models ###################################
################################ Adults ########################################
################################################################################
# Calculate Naive estimate of occupancy
Naive.SS<-sum(SSdata.Adult[3:5], na.rm = T)
Naive.SS/length(SSdata.Adult$Field.Number) #0.4186047

#Split into two years
ADSS2011<-SSdata.Adult[SSdata.Adult$Year=="2011",]
ADSS2016<-SSdata.Adult[SSdata.Adult$Year=="2016",]

#Calculate Naive estimates of occupancy seperately
Naive.SS2011<-sum(ADSS2011[3:5], na.rm = T)
Naive.SS2011/length(ADSS2011$Field.Number) #0.5

#2016
Naive.SS2016<-sum(ADSS2016[3:5], na.rm = T)
Naive.SS2016/length(ADSS2016$Field.Number) #0.3333333

################################################################################
# Making the response variable data frame for the adult occupancy model
################################################################################
SS.frame2<-cbind.data.frame(Haul.1 = SSdata.Adult$Haul.1, 
                            Haul.2 = SSdata.Adult$Haul.2, 
                            Haul.3 = SSdata.Adult$Haul.3)

# Occupancy model data frame using 'unmarked'
SS.umf2<-unmarkedFrameOccu(y=SS.frame2, siteCovs = Habitat)
head(SS.umf2)

################################################################################
# Plot observations
################################################################################
#pdf('SS.juv.detections.pdf',  width=2, height=3)
plot(SS.umf2, panels=1)
#dev.off()

################################################################################
# Model building
################################################################################
fm<-list()
fm[[1]] <- occu(~1 ~1, SS.umf2) #intercept model
#fm[[2]] <- occu(~WV ~1, SS.umf2) # inadequate estimates [Standard errors]
#fm[[3]] <- occu(~Depth ~1, SS.umf2)# inadequate estimates [Standard errors]
#fm[[4]] <- occu(~Submerged ~1, SS.umf2)# inadequate estimates [Standard errors]

# Again we move on with no detection probability covariate
fm[[2]] <- occu(~1 ~Depth, SS.umf2)
fm[[3]] <- occu(~1 ~WV, SS.umf2)
fm[[4]] <- occu(~1 ~Cobble, SS.umf2)

# Model Selection
Mods <- paste("fm", 1:4, sep = "")
Z<-as.data.frame(aictab(cand.set=fm, modnames=Mods, sort=T, second.ord=T, LL=T))
Z

summary(fm[[2]]) 

coef(fm[[2]])
coef(fm[[3]])
coef(fm[[1]])
coef(fm[[4]])

fm2 <- occu(~ 1 ~Depth, SS.umf2) 
(pb <- parboot(fm2, chisq, nsim=999, report=1))

################################################################################
#predictions
Predictpsi<-predict(fm2, type="state", Habitat)
Predictp<-predict(fm2, type="det", Habitat)
head(Predictp)

#Averages
mean(Predictp$Predicted) #0.8007257
sd(Predictp$Predicted)/sqrt(length(Predictp$Predicted)) #0
mean(Predictpsi$Predicted) #0.4223231
sd(Predictpsi$Predicted)/sqrt(length(Predictpsi$Predicted)) #0.03812162

# probability of detecting adult silver shiners at least once after 3 surveys
# given that the site is occupied
1-(1-0.8007257)^3

#################################################################################
# Covariate plot
Depth<-seq(from=min(Habitat$Depth,na.rm=TRUE),
        to=max(Habitat$Depth,na.rm=TRUE), by=0.01)
newdata<-as.data.frame(Depth)
occu<-predict(fm2, type='state', newdata=newdata)

occ.plt2<-ggplot(occu)+
  geom_ribbon(aes(x=Depth,ymin=lower, ymax=upper), fill="grey75")+
  geom_line(aes(x=Depth, y=Predicted), size=0.5)+xlab("Depth (m)")+
  ylim(0,1)+ ylab("Occupancy Probability")+
  ggtitle("Adults")+
  scale_x_continuous(limits=c(0,1.3), # next line ensures 2 decimals on x axis
                     labels=function(x){sprintf("%.2f", x)})+
  theme_me
occ.plt2

################################################################################
########################### Occupancy models ###################################
############################### Combined #######################################
################################################################################
# Calculate Naive estimate of occupancy
Naive.SS<-sum(SSdata.Comb[3:5], na.rm = T)
Naive.SS/length(SSdata.Comb$Field.Number) #0.6976744

#Split into two years
ADSS2011<-SSdata.Comb[SSdata.Comb$Year=="2011",]
ADSS2016<-SSdata.Comb[SSdata.Comb$Year=="2016",]

#Calculate Naive estimates of occupancy separately
Naive.SS2011<-sum(ADSS2011[3:5], na.rm = T)
Naive.SS2011/length(ADSS2011$Field.Number) #0.5

#2016
Naive.SS2016<-sum(ADSS2016[3:5], na.rm = T)
Naive.SS2016/length(ADSS2016$Field.Number) #0.9047619

################################################################################
# Making the response variable data frame for the combined occupancy model
################################################################################
SS.frame3<-cbind.data.frame(Haul.1 = SSdata.Comb$Haul.1, 
                            Haul.2 = SSdata.Comb$Haul.2, 
                            Haul.3 = SSdata.Comb$Haul.3)

# Occupancy model data frame using 'unmarked'
SS.umf3<-unmarkedFrameOccu(y=SS.frame3, siteCovs = Habitat)
head(SS.umf3)

################################################################################
# Plot observations
################################################################################
#pdf('SS.juv.detections.pdf',  width=2, height=3)
plot(SS.umf3, panels=1)
#dev.off()
3/43

################################################################################
# Model building
################################################################################
fm<-list()
fm[[1]] <- occu(~1 ~1, SS.umf3) #intercept model
#fm[[2]] <- occu(~WV ~1, SS.umf3) # inadequate estimates [Standard errors]
fm[[2]] <- occu(~Depth ~1, SS.umf3)# 
fm[[3]] <- occu(~Submerged ~1, SS.umf3)# 

Mods <- paste("fm", 1:3, sep = "")
Z<-as.data.frame(aictab(cand.set=fm, modnames=Mods, sort=T, second.ord=T, LL=T))
Z

# Again we move on with no detection probability covariate
fm[[4]] <- occu(~1 ~Depth, SS.umf3)
fm[[5]] <- occu(~1 ~WV, SS.umf3)
fm[[6]] <- occu(~1 ~Cobble, SS.umf3)
fm[[7]] <- occu(~1 ~Submerged, SS.umf3)

# Model Selection
Mods <- paste("fm", 1:7, sep = "")
Z<-as.data.frame(aictab(cand.set=fm, modnames=Mods, sort=T, second.ord=T, LL=T))
Z

summary(fm[[5]]) 

coef(fm[[5]])
coef(fm[[6]])
coef(fm[[4]])
coef(fm[[1]])
coef(fm[[3]])
coef(fm[[7]])
coef(fm[[2]])

fm5 <- occu(~ 1 ~WV, SS.umf3) 
fm6 <- occu(~ 1 ~Cobble, SS.umf3) 
(pb <- parboot(fm5, chisq, nsim=999, report=1))
(pb <- parboot(fm6, chisq, nsim=999, report=1))

################################################################################
#predictions
Predictpsi<-predict(fm5, type="state", Habitat)
Predictp<-predict(fm5, type="det", Habitat)
head(Predictp)

#Averages
mean(Predictp$Predicted) #0.848892
sd(Predictp$Predicted)/sqrt(length(Predictp$Predicted)) #0
mean(Predictpsi$Predicted) #0.7003215
sd(Predictpsi$Predicted)/sqrt(length(Predictpsi$Predicted)) #0.02588038

# probability of detecting adult silver shiners at least once after 3 surveys
# given that the site is occupied
1-(1-0.848892)^3

################################################################################
#predictions
Predictpsi<-predict(fm6, type="state", Habitat)
Predictp<-predict(fm6, type="det", Habitat)
head(Predictp)

#Averages
mean(Predictp$Predicted) #0.8496829
sd(Predictp$Predicted)/sqrt(length(Predictp$Predicted)) #0
mean(Predictpsi$Predicted) #0.7000564
sd(Predictpsi$Predicted)/sqrt(length(Predictpsi$Predicted)) #0.0242838

# probability of detecting adult silver shiners at least once after 3 surveys
# given that the site is occupied
1-(1-0.8496829)^3

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
  ylim(0,1)+ ylab("Occupancy Probability")+
  ggtitle("Combined")+
  theme_me
occ.plt3

# Covariate plot
Cobble<-seq(from=min(Habitat$Cobble,na.rm=TRUE),
        to=max(Habitat$Cobble,na.rm=TRUE), by=0.01)
newdata<-as.data.frame(Cobble)
occu<-predict(fm6, type='state', newdata=newdata)

occ.plt4<-ggplot(occu)+
  geom_ribbon(aes(x=Cobble,ymin=lower, ymax=upper), fill="grey75")+
  geom_line(aes(x=Cobble, y=Predicted), size=0.5)+
  xlab("Proportion of cobble substrate")+
  ylim(0,1)+ ylab("Occupancy Probability")+
  ggtitle("Combined")+
  theme_me
occ.plt4

#################################################################################
# plot the occupancy covariate plots together
#################################################################################
Figure3<-(occ.plt + occ.plt2)/(occ.plt3+occ.plt4)
 
# export the plots as tiff or png
#tiff("Figure3.tiff", width = 6, height = 6, units = 'in', res = 1000)
#png("Figure3.png", width = 6, height = 6, units = 'in', res = 500)
Figure3 # patchwork makes this really easy to do!
#dev.off()

################################################################################
################################################################################
############################## Community Analysis ##############################
################################################################################
################################################################################
# Make Data frames
Fish.counts2<-as.data.frame(t(aggregate(full.comm[3:32],
                                     by=list(Habitat$Year),sum)))
head(Fish.counts2)
Fish.counts2<-as.data.frame(Fish.counts2[-1,]) #remove top row
colnames(Fish.counts2)<-c("2011","2016") #revise column names

# counts for ggplot
Fish.countsgg<-cbind.data.frame(Species = c(rownames(Fish.counts2),
                                            rownames(Fish.counts2)),
                                   Count = c(Fish.counts2$`2011`, 
                                             Fish.counts2$`2016`),
                                   Year = c(rep("2011",30),
                                            rep("2016",30)))
head(Fish.countsgg)

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
Fish.countsgg.red<-Fish.countsgg[c(1,5,6,8,10,14,15,16,19,21,30,25,
                                   28,29,31,35,36,38,40,44,45,46,
                                   49,51,60,55,58,59),]

# Plot
g1<-ggplot(Fish.countsgg.red, aes(reorder(Species, Relative), 
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
#tiff("Figure3.tiff", width = 6, height = 3, units = 'in', res = 1000)
g1
#dev.off()

################################################################################
################################################################################
## Redundancy Analysis
################################################################################
################################################################################
# transform data
head(ad.SS.comm)

# relative abundance per row
ad.SS.comm.trans<-decostand(ad.SS.comm[3:32], method = "total")
juv.SS.comm.trans<-decostand(juv.SS.comm[3:32], method = "total")
ad.juve.comm.trans<-decostand(ad.juve.comm[3:33], method = "total")
total.comm.trans<-decostand(full.comm[3:32], method = "total")

#choose variables for RDA
head(Habitat)
Habitat.RDA<-cbind.data.frame(Habitat[c(6,7,12,13,15,18,19)])
str(Habitat.RDA)

# Calculate correlations
cor_mat <- cor(Habitat.RDA, method='pearson')

# Generate and export correlation plot
#tiff("Corrplot.tiff", width = 7, height = 7, units = 'in', res = 1000)
corrplot(cor_mat,  method = "ellipse", type = "upper",order = "alphabet", 
         tl.pos='d', tl.cex=0.8, addCoefasPercent=F, tl.col="black")
corrplot(cor_mat, add = TRUE, type = "lower", method = "number", diag = FALSE, 
         order = "alphabet",tl.pos = "n", cl.cex = 0.8,  cl.pos = "n", col="black")
#dev.off()

################################################################################
################################################################################
# Stepwise model selection for juveniles
################################################################################
################################################################################
set.seed(4256)
mod0<-rda(ad.juve.comm.trans ~ 1, Habitat.RDA[1:6]) # intercept
mod1<-rda(ad.juve.comm.trans ~ ., Habitat.RDA[1:6])
rda_select.r <-ordistep(mod0, scope = formula(mod1), direction = "both", 
                        Pin = 0.15, Pout = 0.15, perm.max = 9999)

# rda.fin final model
rda.fin <-rda(ad.juve.comm.trans ~ WV + Cobble + Depth + Submerged, 
              data = Habitat.RDA)
summary(rda.fin)

vif.cca(rda.fin) #variable inflation factor looks good

# Calculate adjusted R2 of RDA
R2adj <- RsquareAdj(rda.fin)$adj.r.squared
R2adj 

# Proportion of variance explained per component
rda.fin$CCA$eig/sum(rda.fin$CCA$eig)

# screeplot
screeplot(rda.fin)

# variable coefficients
coef(rda.fin)

################################################################################
# Triplot
################################################################################
# site scores
scores <- data.frame(Habitat.RDA,rda.fin$CCA$u)

# Species scores
vscores <- data.frame(rda.fin$CCA$v)

# Covariate scores
var.score <- data.frame(rda.fin$CCA$biplot[1:4,1:4]) 
var.score$variables <- c("WV","Cobble","Depth","Submerged")

# Only plotting the most abundant species for visualization purposes
rownames(vscores)
colSums(ad.juve.comm[3:33])
vscores2 <- vscores[c(1,5,6,8,10,14,15,16,19,21,25,28,29,30,31),]
rownames(vscores2)<-c("RoBa","WhSu","BrSt","RaDa", "JoDa", "StSh","CoSh",
                      "SMB","RiCh","RoSh","FaMi","LoDa","CrCh","SSa","SSj")

# Create plot of model
ggRDA <- ggplot(scores, aes(x = RDA1, y = RDA2)) +
  geom_hline(yintercept=0,linetype="dashed",col="black")+
  geom_vline(xintercept=0,linetype="dashed",col="black")+
  geom_point() +
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
                  max.overlaps = 20)+ coord_fixed()+
  geom_text_repel(data = var.score, 
                  aes(x = RDA1, y = RDA2, label = variables),
                  nudge_y = ifelse(var.score$RDA2 > 0, .05, -.05),
                  segment.alpha = 0) +
  labs(x = "RDA Axis 1 (56.0%)", y = "RDA Axis 2 (22.2%)") + 
  theme_me 
  
ggRDA

################################################################################
set.seed(4256)
mod0.1<-rda(total.comm.trans ~ 1, Habitat.RDA[1:6])
mod1.1<-rda(total.comm.trans ~ ., Habitat.RDA[1:6])
rda_select.r.1 <-ordistep(mod0.1, scope = formula(mod1.1), direction = "both", 
                        Pin = 0.15, Pout = 0.15, perm.max = 9999)

# rda.fin final model
rda.fin.1 <-rda(total.comm.trans ~ Depth + Cobble + WV + Submerged, 
                data = Habitat.RDA)
summary(rda.fin.1)

vif.cca(rda.fin.1) #variable inflation factor looks good

# Calculate adjusted R2 of RDA
R2adj.1 <- RsquareAdj(rda.fin.1)$adj.r.squared
R2adj.1 

# Proportion of variance explained per component
rda.fin.1$CCA$eig/sum(rda.fin.1$CCA$eig)

# screeplot
screeplot(rda.fin.1)

# variable coefficients
coef(rda.fin.1)

################################################################################
# Triplot
################################################################################
# site scores
scores.1 <- data.frame(Habitat.RDA,rda.fin.1$CCA$u)

# Species scores
vscores.1 <- data.frame(rda.fin.1$CCA$v)

# Covariate scores
var.score.1 <- data.frame(rda.fin.1$CCA$biplot[1:4,1:4]) 
var.score.1$variables <- c("Depth","Cobble","WV","Submerged")

# Only plotting the most abundant species for visualization purposes
rownames(vscores.1)
colSums(full.comm[3:32])
vscores2.1 <- vscores.1[c(1,5,6,8,10,14,15,16,19,21,25,28,29,30),]
rownames(vscores2.1)<-c("RoBa","WhSu","BrSt","RaDa", "JoDa", "StSh","CoSh",
                      "SMB","RiCh","RoSh","FaMi","LoDa","CrCh","SS")

# Create plot of model
ggRDA2 <- ggplot(scores.1, aes(x = RDA1, y = RDA2)) +
  geom_hline(yintercept=0,linetype="dashed",col="black")+
  geom_vline(xintercept=0,linetype="dashed",col="black")+
  geom_point() +
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
                      label = rownames(vscores2.1)), max.overlaps = 20)+
  geom_text_repel(data = var.score.1, 
                  aes(x = RDA1, y = RDA2, label = variables),
                  nudge_y = ifelse(var.score.1$RDA2 > 0, .05, -.05),
                  segment.alpha = 0) + coord_fixed()+
  labs(x = "RDA Axis 1 (68.9%)", y = "RDA Axis 2 (17.9%)") + 
  theme_me 

ggRDA2

# Export the figure
#png("Figure5.png", width = 5, height = 8, units = 'in', res = 500)
#tiff("Figure5.tiff", width = 5, height = 8, units = 'in', res = 1000)
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

Hauls1  <- Hauls1[3:34] #only include species
Hauls12 <- Hauls12[3:34]
Hauls3  <- Hauls3[3:34]

# Haul 1 remove columns with zero individuals caught
i <- (colSums(Hauls1, na.rm=T) != 0)
Hauls1 <- Hauls1[, i] # keep all the non-zero columns

# Haul2
i <- (colSums(Hauls12, na.rm=T) != 0)
Hauls12 <- Hauls12[, i] 

#Haul3
i <- (colSums(Hauls3, na.rm=T) != 0)
Hauls3 <- Hauls3[, i]

# Remove juvenile and adult differentiation columns
Hauls1  <-Hauls1[-c(25:26)]
Hauls12 <-Hauls12[-c(30:31)]
Hauls3  <-Hauls3[-c(30:31)]

# Transform data for performing PCoA (relative abundance)
Hauls1.trans  <-decostand(Hauls1, method = "total")
Hauls12.trans <-decostand(Hauls12, method = "total")
Hauls3.trans  <-decostand(Hauls3, method = "total")

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

#Test for significant correlations 
set.seed(4256)
protest(X = Hauls12.PCA, Y = Hauls1.PCA, scores="sites", permutations=9999)
protest(X = Hauls3.PCA, Y = Hauls12.PCA, scores="sites", permutations=9999)
protest(X = Hauls3.PCA, Y = Hauls1.PCA, scores="sites", permutations=9999)

# check out residuals for final procrustes
Haul3pro <- procrustes(X = Hauls3.PCA, Y = Hauls1.PCA, 
                       symmetric = TRUE, scale = FALSE)
residuals(Haul3pro)

# procrustes plot
Haul3pro1.df<-cbind.data.frame(Haul3pro$X)
Haul3pro2.df<-cbind.data.frame(Haul3pro$Yrot)
colnames(Haul3pro2.df)<-c( "PC1","PC2")
Haul3progg.df<-rbind.data.frame(Haul3pro1.df,Haul3pro2.df)
Haul3progg.df<-cbind(Haul3progg.df, Site=rep(seq(1:43),2), 
                     ordination = c(rep("Ord3",43),rep("Ord1",43)),
                     Outliers=c(rep("black",10),"red",
                                rep("black",17),"red",
                                rep("black",5),"red","red",
                                "black","red",
                                rep("black",15),"red",rep("black",17),"red",
                                rep("black",5),"red","red",
                                "black","red",rep("black",5)))


#ggplot
# Include sites 36, 29, 35, 11, 38
Haul3progg.df[Haul3progg.df$Site=="36",]
Haul3progg.df[Haul3progg.df$Site=="29",]
Haul3progg.df[Haul3progg.df$Site=="35",]
Haul3progg.df[Haul3progg.df$Site=="11",]
Haul3progg.df[Haul3progg.df$Site=="38",]

Procrustesgg<-ggplot(Haul3progg.df)+
  geom_hline(yintercept=0,linetype="dashed",col="black")+
  geom_vline(xintercept=0,linetype="dashed",col="black")+
  geom_line(aes(x=PC1, y=PC2, group = Site, col=Outliers))+
  scale_colour_manual(values=c("black","red"))+
  geom_point(aes(x=PC1, y=PC2, fill=ordination), shape=21)+
  scale_fill_manual(values=c("black","gray"))+
  theme_me+ 
  annotate("text",x=0.022,y=0.17,label="36")+
  annotate("text",x=0.18,y=-0.08,label="29")+
  annotate("text",x=0.075,y=-0.12,label="35")+
  annotate("text",x=-0.067,y=0.015,label="11")+
  annotate("text",x=0.143,y=-0.07,label="38")+
  xlab("Dimension 1")+ylab("Dimension 2")+
  theme(legend.position = "none")
Procrustesgg

# residuals plot
Residuals.proc<-as.data.frame(cbind(SS=residuals(Haul3pro),
                                    Site=seq(1:43),
                                    Sitecol=c(rep("black",10),"red",
                                              rep("black",17),"red",
                                              rep("black",5),"red","red",
                                              "black","red",
                                              rep("black",5))))

head(Residuals.proc)
str(Residuals.proc)
Residuals.proc$SS<-as.character(Residuals.proc$SS)
Residuals.proc$SS<-as.numeric(Residuals.proc$SS)
Residuals.proc$Site<-as.character(Residuals.proc$Site)
Residuals.proc$Site<-as.numeric(Residuals.proc$Site)
quantile(residuals(Haul3pro))

# ggplot
Residplotgg<-ggplot(Residuals.proc)+
  geom_hline(yintercept=0.04556637  ,linetype="dashed",col="black")+
  geom_hline(yintercept=0.06983697  ,col="black")+
  geom_hline(yintercept=0.11386099  ,linetype="dashed",col="black")+
  geom_col(aes(y=SS, x=Site, fill=Sitecol), width = 0.5)+
  scale_fill_manual(values=c("black","red"))+
  scale_y_continuous(limits=c(0,0.2), labels=function(x){sprintf("%.2f", x)})+
  ylab("Residual deviance")+
  coord_flip()+
  theme_me+
  theme(legend.position = "none")
Residplotgg

# Plot the ordination plot and residual plot together
Figure4<-Procrustesgg+Residplotgg+
  plot_annotation(tag_levels = 'A')+
  plot_layout(widths = c(2, 1))
Figure4 # I think that looks pretty slick

# Export
#png("Figure4.png", width = 7.5, height = 4, units = 'in', res = 500)
#tiff("Figure4.tiff", width = 7.5, height = 4, units = 'in', res = 1000)
Figure4
#dev.off()
