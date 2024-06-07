#### Library ####
library(rgbif) # Get data from GBIF
library(maptools) # For making maps
library(rgeos) # Good for topology
library(raster) # Hadles raster files
library(rgdal) # Helps with geospatial data
library(dismo) # Makes SDMs
library(ggplot2) # Makes nice graphs
library(ggspatial) # Makes nice maps
library(tidyverse) # Nice script - hardly use
library(FlickrAPI) # Communicates with Flickr
library(tiff) # Hadles .tiff files
library(caTools) # For train/test data
library(visreg) # For making GLM plots
library(gridExtra) # For merging polygons together
library(plyr) # For counting duplicates in dataframes
library(chron) # time data
library("rnaturalearth") # For world data for ggspatial
library("rnaturalearthdata") # As above
library(glmmTMB) # For GLMMs
library(MuMIn) # For model analysis
library(corrplot) # For correlation plot
library(psych) # For parallel analyses for PCA
library(vegan) # For PCA stuff
library(MASS) # More PCA stuff
library(ggeffects) # Plotting lines
library(car) # Model analysis
library(bbmle) # Model analysis
library(lme4) # Mixed models
library(multcomp) # Tukey tests
library(dplyr) # Dataframe handling
library(cowplot) # Arranging plots/figures
library(ggsignif) # gg plot stats
library(ggcorrplot) # For corr plots
library(see); library(patchwork) # required for performance
library(performance) # nice for checking model asummptions
library(ggsignif) # For sig box plots
library(ecospat) # Boyce Index

moth <- read.csv("moth_final.csv") ## from figshare



#### Showing how distribution has changed over time ####
## Build map with only nations of interest
data("wrld_simpl") #Get border data from rgeos
myn <- c("United Kingdom", "France", "Germany", "Austria", "Germany", "Denmark",
         "Ireland", "Switzerland", "Belgium", "Netherlands", "Luxembourg", "Czech Republic") ## add or remove "Italy" for sensitivity analyses
my_map <- wrld_simpl[wrld_simpl$NAME %in% myn,]
plot(my_map, axes=T)

## And for ggspatial

world <- ne_countries(scale = "medium", returnclass = "sf")
myn2 <- c("United Kingdom", "France", "Germany", "Austria", "Germany", "Denmark",
          "Ireland", "Switzerland", "Belgium", "Netherlands", "Luxembourg", "Czech Rep.") ## add or remove "Italy" for sensitivity analyses
my_mapg <- world[world$brk_name %in% myn2,]

# And only keep points that are within my_map

motht <- moth

coordinates(motht) <- ~Long+Lat
motht@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
my_map@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

moth$over <- over(motht, my_map)
moth <- moth[complete.cases(moth$over),]

keepov <- c("ID", "year", "Lat", "Long", "source")
moth <- moth[keepov]

# 00s:
# Sample size:
moth00s <- moth[moth$year>1999 & moth$year<2010,]
nrow(moth00s) 

#Map
plot(wrld_simpl, ylim=c(35, 60), xlim=c(-10, 35), axes=T,
     main="00s (n = 830)")
points(moth$Long[moth$year>1999 & moth$year<2010],
       moth$Lat[moth$year>1999 & moth$year<2010], col="red", pch = ".")

# 10s:
# Sample size:
moth10s <- moth[moth$year>2009 & moth$year<2020,]
nrow(moth10s)

#Map
plot(wrld_simpl, ylim=c(35, 60), xlim=c(-10, 35), axes=T,
     main="10s (n = 4565)")
points(moth$Long[moth$year>2009 & moth$year<2020],
       moth$Lat[moth$year>2009 & moth$year<2020], col="blue", pch = ".")

# Also map for Flickr data:
plot(wrld_simpl, ylim=c(35, 60), xlim=c(-10, 35), axes=T,
     main="Flickr (n = 126)")
points(moth$Long[moth$source=="Flickr"],
       moth$Lat[moth$source=="Flickr"], col="red")

# Also map for Twitter data:
plot(wrld_simpl, ylim=c(35, 60), xlim=c(-10, 35), axes=T,
     main="Twitter (n = 21)")
points(moth$Long[moth$source=="Twitter"],
       moth$Lat[moth$source=="Twitter"], col="brown")

# Also map for Instagram data:
plot(wrld_simpl, ylim=c(35, 60), xlim=c(-10, 35), axes=T,
     main="Instagram (n= 147)")
points(moth$Long[moth$source=="Instagram"],
       moth$Lat[moth$source=="Instagram"], col="turquoise")

# Also map for iNaturalist data:
plot(wrld_simpl, ylim=c(35, 60), xlim=c(-10, 35), axes=T,
     main="iNaturalist")
points(moth$Long[moth$source=="iNaturalist"],
       moth$Lat[moth$source=="iNaturalist"], col="dark green")
points(moth$Long[moth$source=="GBIF"],
       moth$Lat[moth$source=="GBIF"], col="red")


#### 00s (Baseline Model) ####

## Climate data edited from worldclim
prec00s <- raster("Climate/00s/prec00s.tif")
covprec0 <- raster("Climate/00s/covprec00s.tif")
tmax00s <- raster("Climate/00s/tmax00s.tif")
covtmax0 <- raster("Climate/00s/covtmax00s.tif")
nl00s <- raster("Climate/00s/nl2012.tif")


### HSM for the 00s ####
hsmmod <- function(x, y) {
  clim <- mask(x, my_map)
  colnames(y) <- c("Long", "Lat")
  obsclim <- raster::extract(x=x, y=y)
  obsclim <- na.omit(obsclim)
  sdm <- bioclim(obsclim)
  pairs(sdm, pa = "p")
  hsm <<- predict(x, sdm)
  plot(hsm, axes = T, ylim=c(30, 60), xlim=c(-10, 40))
  set.seed(sample(0:10000000)[1])
  pa <- randomPoints(x, nrow(y))
  colnames(pa) <- colnames(y)
  absvals <- raster::extract(x, pa)
  print(evaluate(obsclim, absvals, sdm))
  testing.group <- 1
  group.presence <- kfold(x = y, k = 5)
  presence.train <- y[group.presence != testing.group, ]
  presence.test <- y[group.presence == testing.group, ]
  group.background <- kfold(x = pa, k = 5)
  background.train <- pa[group.presence != testing.group, ]
  background.test <- pa[group.presence == testing.group, ]
  bc.model <- bioclim(x = x, p = presence.train)
  predict.presence <<- dismo::predict(object = bc.model, x = x)
  bc.eval <- evaluate(presence.test, background.test, bc.model, x)
  bc.threshold <<- threshold(x = bc.eval, stat = "sensitivity",
                                sensitivity=0.9)
  print(bc.eval@auc)
  
  res <- list(hsm, bc.eval@auc)
  return(res)

}

occur00s <- as.data.frame(cbind(moth$Long[moth$year > 2000 & moth$year < 2009 & moth$source=="GBIF"],
                                 moth$Lat[moth$year > 2000 & moth$year < 2009 & moth$source=="GBIF"]))


hsmodplot <- function(x, y) {
  colnames(y) <- c("Long", "Lat")
  colnames(x) <- c("Long", "Lat")
  ggplot(data = my_mapg) +
    geom_sf(color = "black", fill = "white") +
    ylim(c(35, 60)) +
    xlim(c(-10, 20)) +
    layer_spatial(predict.presence > bc.threshold) +
    scale_fill_binned(low = "white", high = "#999933") +
    xlab("Longitude") +
    ylab("Latitude") +
    geom_sf(color = "black", fill = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(color=c(1,0,1,0)),
          axis.text.y = element_text(color=c(1,0,1,0))) +
    geom_point(data = x, mapping = aes(x = Long, y = Lat), color = "black", size = 1) +
    geom_point(data = y, mapping = aes(x = Long, y = Lat), color = "red", size = 1) +
    theme(legend.position = "none")
}
 ## 10s models ####

prec2010 <- raster("Climate/2010/prec2010.tif")
covprec2010 <- raster("Climate/2010/covprec2010.tif")
tmax2010 <- raster("Climate/2010/tmax2010.tif")
covtmax2010 <- raster("Climate/2010/covtmax2010.tif")
nl2010 <- raster("Climate/2010/nl2012.tif")

clim2010 <- stack(prec2010, covprec2010, tmax2010, covtmax2010, nl2010)
clim2010 <- mask(clim2010, my_map)
occur2010 <- as.data.frame(cbind(moth$Long[moth$year == 2010 & moth$source=="GBIF"],
                                moth$Lat[moth$year == 2010 & moth$source=="GBIF"]))

hsm2010 <- hsmmod(clim2010, occur00s)
hsm2010 <- hsm2010[[1]]
plot(hsm2010)

hsmauc2010 <- replicate(n = 5, hsmmod(clim2010, occur00s), simplify = F)
hsmauc2010m <- mean(c(hsmauc2010[[1]], hsmauc2010[[2]], hsmauc2010[[3]], hsmauc2010[[4]], hsmauc2010[[5]]))
hsmauc2010m

ioccur2010 <- as.data.frame(cbind(moth$Long[moth$year == 2010 & moth$source!="GBIF"],
                                 moth$Lat[moth$year == 2010 & moth$source!="GBIF"]))

hsmodplot(occur2010, ioccur2010)

prec2011 <- raster("Climate/2011/prec2011.tif")
covprec2011 <- raster("Climate/2011/covprec2011.tif")
tmax2011 <- raster("Climate/2011/tmax2011.tif")
covtmax2011 <- raster("Climate/2011/covtmax2011.tif")
nl2011 <- raster("Climate/2011/nl2012.tif")

clim2011 <- stack(prec2011, covprec2011, tmax2011, covtmax2011, nl2011)
clim2011 <- mask(clim2011, my_map)
occur2011 <- as.data.frame(cbind(moth$Long[moth$year == 2011 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2011 & moth$source=="GBIF"]))

hsmauc2011 <- replicate(n = 5, hsmmod(clim2011, occur00s), simplify = F)
hsmauc2011m <- mean(c(hsmauc2011[[1]], hsmauc2011[[2]], hsmauc2011[[3]], hsmauc2011[[4]], hsmauc2011[[5]]))
hsmauc2011m

ioccur2011 <- as.data.frame(cbind(moth$Long[moth$year == 2011 & moth$source!="GBIF"],
                                  moth$Lat[moth$year == 2011 & moth$source!="GBIF"]))

hsmodplot(occur2011, ioccur2011)


prec2012 <- raster("Climate/2012/prec2012.tif")
covprec2012 <- raster("Climate/2012/covprec2012.tif")
tmax2012 <- raster("Climate/2012/tmax2012.tif")
covtmax2012 <- raster("Climate/2012/covtmax2012.tif")
nl2012 <- raster("Climate/2012/nl2012.tif")

clim2012 <- stack(prec2012, covprec2012, tmax2012, covtmax2012, nl2012)
clim2012 <- mask(clim2012, my_map)
occur2012 <- as.data.frame(cbind(moth$Long[moth$year == 2012 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2012 & moth$source=="GBIF"]))

hsmauc2012 <- replicate(n = 5, hsmmod(clim2012, occur00s), simplify = F)
hsmauc2012m <- mean(c(hsmauc2012[[1]], hsmauc2012[[2]], hsmauc2012[[3]], hsmauc2012[[4]], hsmauc2012[[5]]))
hsmauc2012m

ioccur2012 <- as.data.frame(cbind(moth$Long[moth$year == 2012 & moth$source!="GBIF"],
                                  moth$Lat[moth$year == 2012 & moth$source!="GBIF"]))

hsmodplot(occur2012, ioccur2012)


prec2013 <- raster("Climate/2013/prec2013.tif")
covprec2013 <- raster("Climate/2013/covprec2013.tif")
tmax2013 <- raster("Climate/2013/tmax2013.tif")
covtmax2013 <- raster("Climate/2013/covtmax2013.tif")
nl2013 <- raster("Climate/2013/nl2013.tif")

clim2013 <- stack(prec2013, covprec2013, tmax2013, covtmax2013, nl2013)
clim2013 <- mask(clim2013, my_map)
occur2013 <- as.data.frame(cbind(moth$Long[moth$year == 2013 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2013 & moth$source=="GBIF"]))

hsmauc2013 <- replicate(n = 5, hsmmod(clim2013, occur00s), simplify = F)
hsmauc2013m <- mean(c(hsmauc2013[[1]], hsmauc2013[[2]], hsmauc2013[[3]], hsmauc2013[[4]], hsmauc2013[[5]]))
hsmauc2013m

ioccur2013 <- as.data.frame(cbind(moth$Long[moth$year == 2013 & moth$source!="GBIF"],
                                  moth$Lat[moth$year == 2013 & moth$source!="GBIF"]))

hsmodplot(occur2013, ioccur2013)

prec2014 <- raster("Climate/2014/prec2014.tif")
covprec2014 <- raster("Climate/2014/covprec2014.tif")
tmax2014 <- raster("Climate/2014/tmax2014.tif")
covtmax2014 <- raster("Climate/2014/covtmax2014.tif")
nl2014 <- raster("Climate/2014/nl2014.tif")

clim2014 <- stack(prec2014, covprec2014, tmax2014, covtmax2014, nl2014)
clim2014 <- mask(clim2014, my_map)
occur2014 <- as.data.frame(cbind(moth$Long[moth$year == 2014 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2014 & moth$source=="GBIF"]))

hsmauc2014 <- replicate(n = 5, hsmmod(clim2014, occur00s), simplify = F)
hsmauc2014m <- mean(c(hsmauc2014[[1]], hsmauc2014[[2]], hsmauc2014[[3]], hsmauc2014[[4]], hsmauc2014[[5]]))
hsmauc2014m

ioccur2014 <- as.data.frame(cbind(moth$Long[moth$year == 2014 & moth$source!="GBIF"],
                                  moth$Lat[moth$year == 2014 & moth$source!="GBIF"]))

hsmodplot(occur2014, ioccur2014)


prec2015 <- raster("Climate/2015/prec2015.tif")
covprec2015 <- raster("Climate/2015/covprec2015.tif")
tmax2015 <- raster("Climate/2015/tmax2015.tif")
covtmax2015 <- raster("Climate/2015/covtmax2015.tif")
nl2015 <- raster("Climate/2015/nl2015.tif")

clim2015 <- stack(prec2015, covprec2015, tmax2015, covtmax2015, nl2015)
clim2015 <- mask(clim2015, my_map)
occur2015 <- as.data.frame(cbind(moth$Long[moth$year == 2015 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2015 & moth$source=="GBIF"]))

hsmauc2015 <- replicate(n = 5, hsmmod(clim2015, occur00s), simplify = F)
hsmauc2015m <- mean(c(hsmauc2015[[1]], hsmauc2015[[2]], hsmauc2015[[3]], hsmauc2015[[4]], hsmauc2015[[5]]))
hsmauc2015m

ioccur2015 <- as.data.frame(cbind(moth$Long[moth$year == 2015 & moth$source!="GBIF"],
                                  moth$Lat[moth$year == 2015 & moth$source!="GBIF"]))

hsmodplot(occur2015, ioccur2015)


prec2016 <- raster("Climate/2016/prec2016.tif")
covprec2016 <- raster("Climate/2016/covprec2016.tif")
tmax2016 <- raster("Climate/2016/tmax2016.tif")
covtmax2016 <- raster("Climate/2016/covtmax2016.tif")
nl2016 <- raster("Climate/2016/nl2016.tif")

clim2016 <- stack(prec2016, covprec2016, tmax2016, covtmax2016, nl2016)
clim2016 <- mask(clim2016, my_map)
occur2016 <- as.data.frame(cbind(moth$Long[moth$year == 2016 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2016 & moth$source=="GBIF"]))

hsmauc2016 <- replicate(n = 5, hsmmod(clim2016, occur00s), simplify = F)
hsmauc2016m <- mean(c(hsmauc2016[[1]], hsmauc2016[[2]], hsmauc2016[[3]], hsmauc2016[[4]], hsmauc2016[[5]]))
hsmauc2016m

ioccur2016 <- as.data.frame(cbind(moth$Long[moth$year == 2016 & moth$source!="GBIF"],
                                  moth$Lat[moth$year == 2016 & moth$source!="GBIF"]))

hsmodplot(occur2016, ioccur2016)
hsm16 <- hsm


prec2017 <- raster("Climate/2017/prec2017.tif")
covprec2017 <- raster("Climate/2017/covprec2017.tif")
tmax2017 <- raster("Climate/2017/tmax2017.tif")
covtmax2017 <- raster("Climate/2017/covtmax2017.tif")
nl2017 <- raster("Climate/2017/nl2017.tif")

clim2017 <- stack(prec2017, covprec2017, tmax2017, covtmax2017, nl2017)
clim2017 <- mask(clim2017, my_map)
occur2017 <- as.data.frame(cbind(moth$Long[moth$year == 2017 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2017 & moth$source=="GBIF"]))

hsmauc2017 <- replicate(n = 5, hsmmod(clim2017, occur00s), simplify = F)
hsmauc2017m <- mean(c(hsmauc2017[[1]], hsmauc2017[[2]], hsmauc2017[[3]], hsmauc2017[[4]], hsmauc2017[[5]]))
hsmauc2017m

ioccur2017 <- as.data.frame(cbind(moth$Long[moth$year == 2017 & moth$source!="GBIF"],
                                  moth$Lat[moth$year == 2017 & moth$source!="GBIF"]))

hsmodplot(occur2017, ioccur2017)
hsm17 <- hsm


prec2018 <- raster("Climate/2018/prec2018.tif")
covprec2018 <- raster("Climate/2018/covprec2018.tif")
tmax2018 <- raster("Climate/2018/tmax2018.tif")
covtmax2018 <- raster("Climate/2018/covtmax2018.tif")
nl2018 <- raster("Climate/2018/nl2018.tif")

clim2018 <- stack(prec2018, covprec2018, tmax2018, covtmax2018, nl2018)
clim2018 <- mask(clim2018, my_map)
occur2018 <- as.data.frame(cbind(moth$Long[moth$year == 2018 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2018 & moth$source=="GBIF"]))

hsmauc2018 <- replicate(n = 5, hsmmod(clim2018, occur00s), simplify = F)
hsmauc2018m <- mean(c(hsmauc2018[[1]], hsmauc2018[[2]], hsmauc2018[[3]], hsmauc2018[[4]], hsmauc2018[[5]]))
hsmauc2018m

ioccur2018 <- as.data.frame(cbind(moth$Long[moth$year== 2018 & moth$source!="GBIF"],
                                  moth$Lat[moth$year == 2018 & moth$source!="GBIF"]))

hsmodplot(occur2018, ioccur2018)
hsm18 <- hsm

## Boxplots ####
## HS
df2018 <- as.data.frame(cbind(moth$Long[moth$year > 2015 & moth$year<2019],
                              moth$Lat[moth$year > 2015 & moth$year<2019]))
colnames(df2018) <- c("Long", "Lat")

hsmbox <- mean(hsm16, hsm17, hsm18)

hsdf <- raster::extract(x = hsmbox, y = df2018)

df20182 <- as.data.frame(cbind(moth$ID[moth$year > 2015 & moth$year<2019],
                               moth$source[moth$year > 2015 & moth$year<2019],
                               moth$year[moth$year > 2015 & moth$year<2019]))
colnames(df20182) <- c("ID", "Source", "year")

hsdf2 <- cbind(df20182, hsdf)

hsdf2 <- hsdf2[complete.cases(hsdf2$hsdf),]

hsdf2$hsdf <- hsdf2$hsdf + 0.000000001
hsdf2 <- hsdf2[hsdf2$Source!="Twitter",]
table(hsdf2$Source)

hsdf2$loghs <- log(hsdf2$hsdf)

hsdf2$sqhs <- sqrt(hsdf2$hsdf)

summary(hsdf2)

lmhs <- lm(sqhs ~ Source, data=hsdf2)
check_model(lmhs)

summary(lmhs)

anovhs <- aov(sqhs ~ Source, data = hsdf2)
summary(anovhs)

TukeyHSD(anovhs)

gglm2 <- ggplot(data = hsdf2, aes(x = Source, y = sqhs, fill = as.factor(year))) +
  geom_boxplot() +
  ylab("Habitat Suitability") + xlab("Data Source") + theme_bw() + ylim(c(0, 1)) +
  scale_fill_manual(values=c("#88CCEE", "#AA4499","#999933"), name = "Year") +
  theme(text = element_text(size = 24)) +
  geom_signif(comparisons = list(c("Instagram", "GBIF")), map_signif_level = T) +
  geom_signif(comparisons = list(c("Flickr", "iNaturalist")), map_signif_level = T) +
  geom_signif(comparisons = list(c("iNaturalist", "Instagram")), map_signif_level = T)

gglm2

## NL
nlbox <- mean(nl2016, nl2017, nl2018)

nldf <- raster::extract(x = nlbox, y = df2018)

nldf2 <- cbind(df20182, nldf)
nldf2 <- nldf2[nldf2$Source!="Twitter",]
table(nldf2$Source)
nldf2$lognl <- log(nldf2$nldf)

lmnl <- lm(lognl ~ Source, data=nldf2)
check_model(lmnl)

anovnl <- aov(lognl ~ Source, data = nldf2)
summary(anovnl)

TukeyHSD(anovnl)

gglm <- ggplot(data = nldf2, aes(x = Source, y = lognl, fill = as.factor(year))) + geom_boxplot() +
  ylab("Log(Night Light)") + xlab("Data Source") + theme_bw() +
  scale_fill_manual(values=c("#88CCEE", "#AA4499","#999933"), name = "Year") + ylim(c(-2, 6.5)) +
  geom_signif(comparisons = list(c("Flickr", "iNaturalist")), map_signif_level = T) +
  geom_signif(comparisons = list(c("GBIF", "Instagram")), map_signif_level = T) +
  theme(text = element_text(size = 24))

gglm

## Max Temp test for iNaturalist

tmaxbox <- mean(tmax2016, tmax2017, tmax2018)

tmaxdf <- raster::extract(x = tmaxbox, y = df2018)

tmaxdf2 <- cbind(df20182, tmaxdf)
tmaxdf2 <- tmaxdf2[tmaxdf2$Source!="Twitter",]
table(tmaxdf2$Source)
tmaxdf2$logtmax <- log(tmaxdf2$tmaxdf)

lmtmax <- lm(tmaxdf ~ Source, data=tmaxdf2)
check_model(lmtmax)

anovtmax <- aov(tmaxdf ~ Source, data = tmaxdf2)
summary(anovtmax)

TukeyHSD(anovtmax)

ggtmax <- ggplot(data = tmaxdf2, aes(x = Source, y = tmaxdf, fill = as.factor(year))) + geom_boxplot() +
  ylab("Average Maximum Temperature") + xlab("Data Source") + theme_bw() +
  scale_fill_manual(values=c("#88CCEE", "#AA4499","#999933"), name = "Year") +
  geom_signif(comparisons = list(c("Flickr", "GBIF")), map_signif_level = T) +
  geom_signif(comparisons = list(c("iNaturalist", "GBIF")), map_signif_level = T) +
  geom_signif(comparisons = list(c("iNaturalist", "Instagram")), map_signif_level = T) +
  theme(text = element_text(size = 24))

ggtmax

covtbox <- mean(covtmax2016, covtmax2017, covtmax2018)

covtdf <- raster::extract(x = covtbox, y = df2018)

covtdf2 <- cbind(df20182, covtdf)
covtdf2 <- covtdf2[covtdf2$Source!="Twitter",]
table(covtdf2$Source)
covtdf2$logcovt <- log(covtdf2$covtdf)

lmcovt <- lm(covtdf ~ Source, data=covtdf2)
check_model(lmcovt)

anovcovt <- aov(covtdf ~ Source, data = covtdf2)
summary(anovcovt)

TukeyHSD(anovcovt)

ggcovt <- ggplot(data = covtdf2, aes(x = Source, y = covtdf, fill = as.factor(year))) + geom_boxplot() +
  ylab("Average Maximum Temperature") + xlab("Data Source") + theme_bw() +
  scale_fill_manual(values=c("#88CCEE", "#AA4499","#999933"), name = "Year") +
  geom_signif(comparisons = list(c("Flickr", "GBIF")), map_signif_level = T) +
  geom_signif(comparisons = list(c("iNaturalist", "GBIF")), map_signif_level = T) +
  geom_signif(comparisons = list(c("iNaturalist", "Instagram")), map_signif_level = T) +
  theme(text = element_text(size = 24))

ggcovt

precbox <- mean(prec2016, prec2017, prec2018)

precdf <- raster::extract(x = precbox, y = df2018)

precdf2 <- cbind(df20182, precdf)
precdf2 <- precdf2[precdf2$Source!="Twitter",]
table(precdf2$Source)
precdf2$logprec <- log(precdf2$precdf)

lmprec <- lm(precdf ~ Source, data=precdf2)
check_model(lmprec)

anovprec <- aov(precdf ~ Source, data = precdf2)
summary(anovprec)

TukeyHSD(anovprec)

ggprec <- ggplot(data = precdf2, aes(x = Source, y = precdf, fill = as.factor(year))) + geom_boxplot() +
  ylab("Average Maximum Temperature") + xlab("Data Source") + theme_bw() +
  scale_fill_manual(values=c("#88CCEE", "#AA4499","#999933"), name = "Year") +
  geom_signif(comparisons = list(c("Flickr", "GBIF")), map_signif_level = T) +
  geom_signif(comparisons = list(c("iNaturalist", "GBIF")), map_signif_level = T) +
  geom_signif(comparisons = list(c("iNaturalist", "Instagram")), map_signif_level = T) +
  theme(text = element_text(size = 24))

ggprec

#### Recorder effort ####

### Reset occur ####

occura2010 <- as.data.frame(cbind(moth$Long[moth$year==2010],
                                  moth$Lat[moth$year==2010]))
occura2011 <- as.data.frame(cbind(moth$Long[moth$year==2011],
                                  moth$Lat[moth$year==2011]))
occura2012 <- as.data.frame(cbind(moth$Long[moth$year==2012],
                                  moth$Lat[moth$year==2012]))
occura2013 <- as.data.frame(cbind(moth$Long[moth$year==2013],
                                  moth$Lat[moth$year==2013]))
occura2014 <- as.data.frame(cbind(moth$Long[moth$year==2014],
                                  moth$Lat[moth$year==2014]))
occura2015 <- as.data.frame(cbind(moth$Long[moth$year==2015],
                                  moth$Lat[moth$year==2015]))
occura2016 <- as.data.frame(cbind(moth$Long[moth$year==2016],
                                  moth$Lat[moth$year==2016]))
occura2017 <- as.data.frame(cbind(moth$Long[moth$year==2017],
                                  moth$Lat[moth$year==2017]))
occura2018 <- as.data.frame(cbind(moth$Long[moth$year==2018],
                                  moth$Lat[moth$year==2018]))


#### Blackbird recorder effort ####
bb <- readOGR("EBBA2_Turdus_merula_2020_06_30/EBBA2_Turdus_merula_2020_06_30.shp")
bb$abundance_2 <- ifelse(bb$abundance_=="A", mean(c(1, 9)), 
                         ifelse(bb$abundance_=="B", mean(c(10, 99)),
                                ifelse(bb$abundance_=="C", mean(c(100, 999)),
                                       ifelse(bb$abundance_=="D", mean(c(1000, 9999)),
                                              ifelse(bb$abundance_=="E", mean(c(10000, 99999)),
                                                     ifelse(bb$abundance_=="F", 100000, bb$abundance_))))))
bb <- spTransform(bb, CRS("+proj=longlat +datum=WGS84"))
bb$abundance_2 <- as.numeric(bb$abundance_2)

bb@proj4string
crop@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

plot(my_map, main = "Blackbird Abundance Test", axes=T)
plot(bb, col = bb$abundance_2, add = T)

#### iNaturalist RE ####
bbi <- read.csv("bbREiNat.csv", header = T, stringsAsFactors = F)
colnames(bbi) [c(1, 3, 22, 23)] <- c("ID", "year", "Lat", "Long")
keep5 <- c("ID", "year", "Lat", "Long")
bbi <- bbi[keep5]
bbi$year <- substr(bbi$year, 0, 4)
bbi$Source <- "iNaturalist"
bbi <- bbi[bbi$year>1999 & bbi$year<2019,]
bbi <- bbi[complete.cases(bbi$Lat),]
bbi <- bbi[complete.cases(bbi$Long),]

#### Flickr RE ####
bbF <- read.csv("bbREFlickr.csv", header = T, stringsAsFactors = F)
colnames(bbF) [c(3, 4, 5, 6)] <- c("ID", "Lat", "Long", "year")
bbF <- bbF[keep5]
bbF$year <- substr(bbF$year, 0, 4)
bbF$Source <- "Flickr"
bbF <- bbF[bbF$year>1999 & bbF$year<2019,]
bbF <- bbF[complete.cases(bbF$Lat),]
bbF <- bbF[complete.cases(bbF$Long),]

#### Instagram RE ####
bbI <- read.csv("bbREInsta.csv", header = T, stringsAsFactors = F)
head(bbI)
bbI <- bbI[keep5]
bbI$source <- "Instagram"

#### Twitter RE ####
bbT <- read.csv("bbRETwit.csv", header = T, stringsAsFactors = F)
head(bbT)
bbT <- bbT[keep5]
bbT$source <- "Twitter"

#### GBIF RE ####

bbproc <- function(x) {
  bb <- read.csv(x, header = T, stringsAsFactors = F, fill = T)
  colnames(bb)[c(1, 133, 134)] <- c("ID", "Lat", "Long")
  bb$inatb <- substr(bb$references, 13, 23)
  bb <- bb[bb$inatb!="inaturalist",]
  keep <- c("ID", "year", "Lat", "Long")
  bb <- bb[keep]
  bb$Lat <- as.numeric(bb$Lat)
  bb$year <- as.numeric(bb$year)
  bb <- bb[complete.cases(bb$Lat),]
  bb <- bb[complete.cases(bb$Long),]
  bb <- bb[complete.cases(bb$year),]
  bb <- bb[bb$year>1999 & bb$year<2019,]
  bb$source <- "GBIF"
  
  bbt <- bb
  
  coordinates(bbt) <- ~Long+Lat
  bbt@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  
  bb$over <- over(bbt, my_map)
  bb <- bb[complete.cases(bb$over),]
  
  keepov <- c("ID", "year", "Lat", "Long", "source")
  bb <- bb[keepov]
}

bb1 <- bbproc("bbREGBIF1.csv")
bb2 <- bbproc("bbREGBIF2.csv")
bb3 <- bbproc("bbREGBIF3.csv")
bb4 <- bbproc("bbREGBIF4.csv")
bb5 <- bbproc("bbREGBIF5.csv")
bb6 <- bbproc("bbREGBIF6.csv")
bb7 <- bbproc("bbREGBIF7.csv")
bb8 <- bbproc("bbREGBIF8.csv")
bb9 <- bbproc("bbREGBIF9.csv")
bb10 <- bbproc("bbREGBIF10.csv")
bb11 <- bbproc("bbREGBIF11.csv")

bbG <- rbind(bb1, bb2, bb3, bb4, bb5, bb6, bb7, bb8, bb9, bb10, bb11)
summary(bbG)

#### Making RE dataframes for models ####

bbrast <- function(x) {
  bbG <- x
  keep6 <- c("Long", "Lat")
  bbG <- bbG[keep6]
  coordinates(bbG) <- ~Long+Lat
  rast <- raster(ncol = ncol(prec2010), nrow = nrow(prec2010))
  extent(rast) <- extent(prec2010)
  bbG <- rasterize(bbG, rast, fun = sum)
  bbG <- mask(bbG, my_map)
  
}

bbG2016 <- bbrast(bbG[bbG$year == 2016,])
bbG2016[is.na(bbG2016)] <- 0

bb@data$RE.G16 <- raster::extract(bbG2016, 
                                  bb,
                                  fun = sum)
bb@data$RE.G16 <- bb@data$RE.G16/bb@data$abundance_2


bbi2016 <- bbrast(bbi[bbi$year == 2016,])
bbi2016[is.na(bbi2016)] <- 0

bb@data$RE.i16 <- raster::extract(bbi2016, 
                                  bb,
                                  fun = sum)
bb@data$RE.i16 <- bb@data$RE.i16/bb@data$abundance_2

bbF$Lat <- as.numeric(bbF$Lat)
bbF$Long <- as.numeric(bbF$Long)
bbF <- bbF[complete.cases(bbF$Lat),]

bbF2016 <- bbrast(bbF[bbF$year == 2016,])
bbF2016[is.na(bbF2016)] <- 0

bb@data$RE.F16 <- raster::extract(bbF2016, 
                                  bb,
                                  fun = sum)
bb@data$RE.F16 <- bb@data$RE.F16/bb@data$abundance_2

bbI2016 <- bbrast(bbI[bbI$year == 2016,])
bbI2016[is.na(bbI2016)] <- 0

bb@data$RE.I16 <- raster::extract(bbI2016, 
                                  bb,
                                  fun = sum)
bb@data$RE.I16 <- bb@data$RE.I16/bb@data$abundance_2

##

REdf16 <- raster::extract(x = bb, y = occura2016)
REdf16a <- raster::extract(x = hsm16, y = occura2016)
REdf16a <- as.data.frame(REdf16a)
colnames(REdf16a) <- "HS"
REdf16 <- cbind(REdf16, REdf16a)
keep9 <- c("HS", "RE.G16", "RE.i16", "RE.F16", "RE.I16")
REdf16 <- REdf16[keep9]
REdf16$RE.G16 <- ifelse(REdf16$RE.G16==Inf, 0, REdf16$RE.G16)
REdf16$RE.i16 <- ifelse(REdf16$RE.i16==Inf, 0, REdf16$RE.i16)
REdf16$RE.F16 <- ifelse(REdf16$RE.F16==Inf, 0, REdf16$RE.F16)
REdf16$RE.I16 <- ifelse(REdf16$RE.I16==Inf, 0, REdf16$RE.I16)

RE2016 <- stack(hsm16, bbG2016, bbi2016, bbF2016, bbI2016)

set.seed(2016)
rp16 <- randomPoints(RE2016, nrow(REdf16))
rp16 <- as.data.frame(rp16)
colnames(rp16) <- colnames(occura2016)
points16 <- rp16
rp16 <- raster::extract(x = bb, y = rp16)
rp16a <- raster::extract(x = hsm16, y = points16)
rp16a <- as.data.frame(rp16a)
colnames(rp16a) <- "HS"
rp16 <- cbind(rp16, rp16a)
rp16 <- rp16[keep9]
rp16$RE.G16 <- ifelse(rp16$RE.G16==Inf, 0, rp16$RE.G16)
rp16$RE.i16 <- ifelse(rp16$RE.i16==Inf, 0, rp16$RE.i16)
rp16$RE.F16 <- ifelse(rp16$RE.F16==Inf, 0, rp16$RE.F16)
rp16$RE.I16 <- ifelse(rp16$RE.I16==Inf, 0, rp16$RE.I16)

sq16 <- rbind(REdf16, rp16)
string16 <- c(rep(1, nrow(REdf16)), rep(0, nrow(rp16)))

REdf16 <- data.frame(cbind(pa=string16, sq16))
colnames(REdf16)

occurb2016 <- as.data.frame(cbind(moth$source[moth$year==2016]))
occurb2016 <- rbind(occurb2016, occurb2016)
REdf16 <- cbind(REdf16, occurb2016)
colnames(REdf16) [7] <- "source"

##

bbG2017 <- bbrast(bbG[bbG$year == 2017,])
bbG2017[is.na(bbG2017)] <- 0

bb@data$RE.G17 <- raster::extract(bbG2017, 
                                  bb,
                                  fun = sum)
bb@data$RE.G17 <- bb@data$RE.G17/bb@data$abundance_2


bbi2017 <- bbrast(bbi[bbi$year == 2017,])
bbi2017[is.na(bbi2017)] <- 0

bb@data$RE.i17 <- raster::extract(bbi2017, 
                                  bb,
                                  fun = sum)
bb@data$RE.i17 <- bb@data$RE.i17/bb@data$abundance_2

bbF2017 <- bbrast(bbF[bbF$year == 2017,])
bbF2017[is.na(bbF2017)] <- 0

bb@data$RE.F17 <- raster::extract(bbF2017, 
                                  bb,
                                  fun = sum)
bb@data$RE.F17 <- bb@data$RE.F17/bb@data$abundance_2

bbI2017 <- bbrast(bbI[bbI$year == 2017,])
bbI2017[is.na(bbI2017)] <- 0

bb@data$RE.I17 <- raster::extract(bbI2017, 
                                  bb,
                                  fun = sum)
bb@data$RE.I17 <- bb@data$RE.I17/bb@data$abundance_2

##

REdf17 <- raster::extract(x = bb, y = occura2017)
REdf17a <- raster::extract(x = hsm17, y = occura2017)
REdf17a <- as.data.frame(REdf17a)
colnames(REdf17a) <- "HS"
REdf17 <- cbind(REdf17, REdf17a)
keep9 <- c("HS", "RE.G17", "RE.i17", "RE.F17", "RE.I17")
REdf17 <- REdf17[keep9]
REdf17$RE.G17 <- ifelse(REdf17$RE.G17==Inf, 0, REdf17$RE.G17)
REdf17$RE.i17 <- ifelse(REdf17$RE.i17==Inf, 0, REdf17$RE.i17)
REdf17$RE.F17 <- ifelse(REdf17$RE.F17==Inf, 0, REdf17$RE.F17)
REdf17$RE.I17 <- ifelse(REdf17$RE.i17==Inf, 0, REdf17$RE.I17)

RE2017 <- stack(hsm17, bbG2017, bbi2017, bbF2017, bbI2017)
nrow(REdf17)

rp17 <- randomPoints(RE2017, nrow(REdf17))
rp17 <- as.data.frame(rp17)
colnames(rp17) <- colnames(occura2017)
points17 <- rp17
rp17 <- raster::extract(x = bb, y = rp17)
rp17a <- raster::extract(x = hsm17, y = points17)
rp17a <- as.data.frame(rp17a)
colnames(rp17a) <- "HS"
rp17 <- cbind(rp17, rp17a)
rp17 <- rp17[keep9]
rp17$RE.G17 <- ifelse(rp17$RE.G17==Inf, 0, rp17$RE.G17)
rp17$RE.i17 <- ifelse(rp17$RE.i17==Inf, 0, rp17$RE.i17)
rp17$RE.F17 <- ifelse(rp17$RE.F17==Inf, 0, rp17$RE.F17)
rp17$RE.I17 <- ifelse(rp17$RE.I17==Inf, 0, rp17$RE.I17)

sq17 <- rbind(REdf17, rp17) 
string17 <- c(rep(1, nrow(REdf17)), rep(0, nrow(rp17)))

REdf17 <- data.frame(cbind(pa=string17, sq17))
colnames(REdf17)

occurb2017 <- as.data.frame(cbind(moth$source[moth$year==2017]))
occurb2017 <- rbind(occurb2017, occurb2017)
REdf17 <- cbind(REdf17, occurb2017)
colnames(REdf17) [7] <- "source"

##

bbG2018 <- bbrast(bbG[bbG$year == 2018,])
bbG2018[is.na(bbG2018)] <- 0

bb@data$RE.G18 <- raster::extract(bbG2018, 
                                  bb,
                                  fun = sum)
bb@data$RE.G18 <- bb@data$RE.G18/bb@data$abundance_2


bbi2018 <- bbrast(bbi[bbi$year == 2018,])
bbi2018[is.na(bbi2018)] <- 0

bb@data$RE.i18 <- raster::extract(bbi2018, 
                                  bb,
                                  fun = sum)
bb@data$RE.i18 <- bb@data$RE.i18/bb@data$abundance_2

bbF2018 <- bbrast(bbF[bbF$year == 2018,])
bbF2018[is.na(bbF2018)] <- 0

bb@data$RE.F18 <- raster::extract(bbF2018, 
                                  bb,
                                  fun = sum)
bb@data$RE.F18 <- bb@data$RE.F18/bb@data$abundance_2

bbI2018 <- bbrast(bbI[bbI$year == 2018,])
bbI2018[is.na(bbI2018)] <- 0

bb@data$RE.I18 <- raster::extract(bbI2018, 
                                  bb,
                                  fun = sum)
bb@data$RE.I18 <- bb@data$RE.I18/bb@data$abundance_2

##

REdf18 <- raster::extract(x = bb, y = occura2018)
REdf18a <- raster::extract(x = hsm18, y = occura2018)
REdf18a <- as.data.frame(REdf18a)
colnames(REdf18a) <- "HS"
REdf18 <- cbind(REdf18, REdf18a)
keep9 <- c("HS", "RE.G18", "RE.i18", "RE.F18", "RE.I18")
REdf18 <- REdf18[keep9]
REdf18$RE.G18 <- ifelse(REdf18$RE.G18==Inf, 0, REdf18$RE.G18)
REdf18$RE.i18 <- ifelse(REdf18$RE.i18==Inf, 0, REdf18$RE.i18)
REdf18$RE.F18 <- ifelse(REdf18$RE.F18==Inf, 0, REdf18$RE.F18)
REdf18$RE.I18 <- ifelse(REdf18$RE.I18==Inf, 0, REdf18$RE.I18)

RE2018 <- stack(hsm18, bbG2018, bbi2018, bbF2018, bbI2018)

set.seed(2018)
rp18 <- randomPoints(RE2018, nrow(REdf18))
rp18 <- as.data.frame(rp18)
colnames(rp18) <- colnames(occura2018)
points18 <- rp18
rp18 <- raster::extract(x = bb, y = rp18)
rp18a <- raster::extract(x = hsm18, y = points18)
rp18a <- as.data.frame(rp18a)
colnames(rp18a) <- "HS"
rp18 <- cbind(rp18, rp18a)
rp18 <- rp18[keep9]
rp18$RE.G18 <- ifelse(rp18$RE.G18==Inf, 0, rp18$RE.G18)
rp18$RE.i18 <- ifelse(rp18$RE.i18==Inf, 0, rp18$RE.i18)
rp18$RE.F18 <- ifelse(rp18$RE.F18==Inf, 0, rp18$RE.F18)
rp18$RE.I18 <- ifelse(rp18$RE.I18==Inf, 0, rp18$RE.I18)

sq18 <- rbind(REdf18, rp18) 
string18 <- c(rep(1, nrow(REdf18)), rep(0, nrow(rp18)))

REdf18 <- data.frame(cbind(pa=string18, sq18))
colnames(REdf18)

occurb2018 <- as.data.frame(cbind(moth$source[moth$year==2018]))
occurb2018 <- rbind(occurb2018, occurb2018)
REdf18 <- cbind(REdf18, occurb2018)
colnames(REdf18) [7] <- "source"

##

colnames(REdf16) [c(3, 4, 5, 6, 7)] <- c("RE.G", "RE.i", "RE.F", "RE.I", "source")
colnames(REdf17) [c(3, 4, 5, 6)] <- c("RE.G", "RE.i", "RE.F", "RE.I")
colnames(REdf18) [c(3, 4, 5, 6, 7)] <- c("RE.G", "RE.i", "RE.F", "RE.I",
                                            "source")

REdf <- rbind(REdf16, REdf17, REdf18)
REdf$RE.G <- ifelse(is.na(REdf$RE.G), 0, REdf$RE.G)
REdf$RE.i <- ifelse(is.na(REdf$RE.i), 0, REdf$RE.i)
REdf$RE.F <- ifelse(is.na(REdf$RE.F), 0, REdf$RE.F)
REdf$RE.I <- ifelse(is.na(REdf$RE.I), 0, REdf$RE.I)

REdfG <- REdf[REdf$source=="GBIF",]
REdfi <- REdf[REdf$source=="iNaturalist",]
REdfF <- REdf[REdf$source=="Flickr",]
REdfI <- REdf[REdf$source=="Instagram",]

REdfG$mRE <- mean(REdfG$RE.G)
REdfG$sdRE <- sd(REdfG$RE.G)
REdfG$sRE <- (REdfG$RE.G-REdfG$mRE)/REdfG$sdRE

REdfi$mRE <- mean(REdfi$RE.i)
REdfi$sdRE <- sd(REdfi$RE.i)
REdfi$sRE <- (REdfi$RE.i-REdfi$mRE)/REdfi$sdRE

REdfF$mRE <- mean(REdfF$RE.F)
REdfF$sdRE <- sd(REdfF$RE.F)
REdfF$sRE <- (REdfF$RE.F-REdfF$mRE)/REdfF$sdRE

REdfI$mRE <- mean(REdfI$RE.I)
REdfI$sdRE <- sd(REdfI$RE.I)
REdfI$sRE <- (REdfI$RE.I-REdfI$mRE)/REdfI$sdRE

REdfG <- REdfG[complete.cases(REdfG$HS),]
REdfi <- REdfi[complete.cases(REdfi$HS),]
REdfF <- REdfF[complete.cases(REdfF$HS),]
REdfI <- REdfI[complete.cases(REdfI$HS),]

options(na.action = "na.fail")

#### Modeling RE ####

## GBIF recorder effort

glmmG <- glmmTMB(pa ~ HS + sRE,
                 family=binomial(link = "logit"), data = REdfG)
dredge(glmmG)

ggHS <- plot(ggpredict(glmmG, c("HS [all]")), add = T,
             dot.size = 0.7, show.title = F)+
  labs(x = "Habitat Suitability", y = "Presence/Absence", col = "Year") +
  scale_y_continuous(breaks = c(0, 1)) + theme_bw() + ggtitle("GBIF")

ggRE <- plot(ggpredict(glmmG, c("sRE [all]")), add = T,
             dot.size = 0.7, show.title = F) +
  labs(x = "Recorder Effort", y = "Presence/Absence", col = "Year") +
  scale_y_continuous(breaks = c(0, 1)) + xlim(c(-1, 10)) + theme_bw()

plot_grid(ggHS, ggRE, nrow = 2, labels = c("(a)", "(b)"),
          label_size = 10)
# GBIF has a positive relationship with HS and RE

## iNaturalist recorder effort
glmmi <- glmmTMB(pa ~ HS + sRE,
                 family=binomial(link = "logit"), data = REdfi)
dredge(glmmi)

ggHSi <- plot(ggpredict(glmmi, c("HS [all]")), add = T,
              dot.size = 0.3, show.title = F) +
  labs(x = "Habitat Suitability", y = "Presence/Absence", col = "Year") +
  scale_y_continuous(breaks = c(0, 1)) + theme_bw() + ggtitle("iNaturalist")

ggREi <- plot(ggpredict(glmmi, c("sRE [all]")), add = T,
              dot.size = 0.7, show.title = F) +
  labs(x = "Recorder Effort", y = "Presence/Absence", col = "Year") +
  scale_y_continuous(breaks = c(0, 1)) + xlim(c(-0.5, 3)) + theme_bw()

plot_grid(ggHSi, ggREi, nrow = 2, labels = c("(a)", "(b)"),
          label_size = 10)
# iNaturalist has a positive relationship with HS and RE

## Flickr recorder effort
glmmF <- glmmTMB(pa ~ HS + sRE,
                 family=binomial(link = "logit"), data = REdfF)
dredge(glmmF)

ggHSF <- plot(ggpredict(glmmF, c("HS [all]")), add = T,
              dot.size = 0.7, show.title = F) +
  labs(x = "Habitat Suitability", y = "Presence/Absence", col = "Year") +
  scale_y_continuous(breaks = c(0, 1)) + theme_bw() + ggtitle("Flickr")

ggREF <- plot(ggpredict(glmmF, c("sRE [all]")), add = T,
              dot.size = 0.7, show.title = F) +
  labs(x = "Recorder Effort", y = "Presence/Absence", col = "Year") +
  scale_y_continuous(breaks = c(0, 1)) + xlim(c(-1, 4)) + theme_bw()

plot_grid(ggHSF, ggREF, nrow = 2, labels = c("(a)", "(b)"),
          label_size = 10)

# Flickr has a positive relationship with HS and RE

## Instagram recorder effort

glmmI <- glmmTMB(pa ~ HS + sRE,
                 family=binomial(link = "logit"), data = REdfI)
dredge(glmmI)

ggHSI <- plot(ggpredict(glmmI, c("HS [all]")), add = T,
              dot.size = 0.7, show.title = F) +
  labs(x = "Habitat Suitability", y = "Presence/Absence", col = "Year") +
  scale_y_continuous(breaks = c(0, 1)) + theme_bw() + ggtitle("Instagram")

ggREI <- plot(ggpredict(glmmI, c("sRE [all]")), add = T,
              dot.size = 0.7, show.title = F) +
  labs(x = "Recorder Effort", y = "Presence/Absence", col = "Year") +
  scale_y_continuous(breaks = c(0, 1)) + xlim(c(-0.5, 1.1)) + theme_bw()

plot_grid(ggHSI, ggREI, nrow = 2, labels = c("(a)", "(b)"),
          label_size = 10)

# Instagram has a negative relationship with HS and no relationship
# with RE

plot_grid(ggHS, ggHSi, ggHSF, ggHSI, 
          ggRE, ggREi, ggREF, ggREI,
          nrow = 2,
          labels = c("(a)", "(c)", "(e)", "(g)",
                     "(b)", "(d)", "(f)", "(h)",
                     label_size = 3), label_x = -0.04,
          greedy=T)
