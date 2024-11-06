############################################
###### packages and data processing ########
## data on figshare - download and ignore ##
############################################


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

#### Importing GBIF data ####
Eup <- occ_data(scientificName = "Euplagia quadripunctaria", hasCoordinate = TRUE, limit = 200000)
Eup <- Eup$data
moth1 <- Eup[ !(Eup$lifeStage %in% "LARVA"), ] # Remove larvae
moth1 <- moth1[!(moth1$individualCount %in% 0), ] # Remove any absence data recorded

### Separating iNaturalist data ## for iNat sensitivity analyses
#moth1$inatb <- substr(moth1$references, 13, 23)
# iNaturalist is included in the references column as a hyperlink
#moth1 <- moth1[moth1$inatb!="inaturalist",]
# remove this to avoid duplicating iNaturalist within GBIF

#imoth <- read.csv("iNatMoth.csv", header = T, stringsAsFactors = F)
# Read in dataframe downloaded from iNaturalist

#### Flickr shortcut ####
## Sourced from Flickr API code, which accesses Flickr via an API and downloads data.
moth2 <- read.csv("FlickMoth.csv", header = T, stringsAsFactors = F)
summary(moth2)

moth2.5 <- read.csv("FlickMoth2.csv", header = T, stringsAsFactors = F) # For other names of the moth
moth2 <- rbind(moth2, moth2.5)

#### Importing manually-gathered Twitter data ####
# Manually checked life stage and date
moth3 <- read.csv("twitter.csv", header=T, stringsAsFactors = F)
summary(moth3)
head(moth3)

#### Importing manually-gathered Instagram data ####
moth4 <- read.csv("insta.csv", header= T, stringsAsFactors = F)
summary(moth4)

#### Processing GBIF data ####
# I only want to keep a few of these
colnames(moth1)[c(1, 3, 4)] <- c("ID", "Lat", "Long")
keep <- c("ID", "year", "Lat", "Long")
moth1 <- moth1[keep]
moth1 <- moth1[complete.cases(moth1$Lat),]
moth1 <- moth1[complete.cases(moth1$Long),]
moth1 <- moth1[complete.cases(moth1$year),]
moth1$source <- "GBIF" # Denote the source to compare with other sources of data

#### Processing iNaturalist data ####

# Unlike with GBIF, iNaturalist does provide the actual date of the photo, so I need to
# extract the data that are within my range
getSeason <- function(DATES) {
  NS <- as.Date("2012-10-1", format = "%Y-%m-%d") # Not Summer
  S <- as.Date("2012-7-1",  format = "%Y-%m-%d") # Summer
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= NS | d < S, "Not Summer",
          "Summer")
}

#colnames(imoth) [c(1, 3, 22, 23)] <- c("ID", "year", "Lat", "Long")
#imoth <- imoth[keep]
#imoth <- imoth[complete.cases(imoth$Lat),]
#imoth <- imoth[complete.cases(imoth$Long),]
#imoth <- imoth[complete.cases(imoth$year),]
#imoth$Season <- getSeason(imoth$year)
#imoth <- imoth[imoth$Season=="Summer",]
#imoth$year <- substr(imoth$year, 0, 4)
#imoth$source <- "iNaturalist"

#### Processing Flickr data ####



# Then addin a colum  with the results
moth2$Season <- getSeason(moth2$date)

# And remove those I don't want
moth2 <- moth2[moth2$Season=="Summer",]

# Then, once again, streamline the data
moth2$year <- substr(moth2$date, 0, 4)
colnames(moth2)[c(1, 4, 5)] <- c("ID", "Lat", "Long")
keep2 <- c("ID", "year", "Lat", "Long")
moth2 <- moth2[keep]
moth2$source <- "Flickr"

#### Processing Twitter data ####
# Similar to above
colnames(moth3)[2] <- "year"
keep3 <- c("ID", "year", "Lat", "Long")
moth3 <- moth3[keep3]
moth3$source <- "Twitter"

#### Processing Instagram data ####
colnames(moth4)
keep4 <- c("ID", "year", "Lat", "Long")
moth4 <- moth4[keep4]
moth4$source <- "Instagram"

#### Merge all into one big, mothy dataframe ####
moth.5 <- merge(moth1, moth2, all=T)
moth.6 <- merge(moth.5, moth3, all=T)
#moth.7 <- merge(moth.6, imoth, all=T)
moth <- merge(moth.6, moth4, all=T)
moth$year <- as.numeric(moth$year)
moth$Lat <- as.numeric(moth$Lat)
moth$Long <- as.numeric(moth$Long)
moth <- moth[ !(moth$ID %in% 1893201994), ] # Remove datapoint from Antarctica
moth <- moth[!duplicated(moth[,c("year", 'Lat', 'Long')]),] # Remove duplicates
moth <- moth[moth$year>1999 & moth$year<2019,]
moth <- moth[complete.cases(moth$Long),]
nrow(moth)

#write.csv(moth, "moth_final.csv")

#### Descriptive statistics ####
# Histograms based on timebin size
qplot(moth$year, geom="histogram", fill=moth$source, col=I("black"), binwidth=1) + 
  xlab("Year") + ylab("No. Occurances") + 
  theme_minimal() + labs(fill = "Data Source") + 
  scale_fill_manual(values = c("#d8b365", "#5ab4ac", "white", "brown", "yellow")) +
  theme_bw()

# For Flickr alone:
nrow(moth[moth$source == "Flickr",]) # 135
qplot(moth$year[moth$source=="Flickr"], geom="histogram", fill=I("#5ab4ac"),
      col=I("white"), binwidth=1) + ggtitle("Flickr (n= 126)") + xlab("Year") +
  ylab("No. Occurances")

# For Twitter alone:
nrow(moth[moth$source == "Twitter",]) #16
qplot(moth$year[moth$source=="Twitter"], geom="histogram", fill=I("#5ab4ac"),
      col=I("white"), binwidth=1) + ggtitle("Twitter (n = 21)") + xlab("Year") +
  ylab("No. Occurances")

# For Instagram alone:
nrow(moth[moth$source == "Instagram",]) # 148
qplot(moth$year[moth$source=="Instagram"], geom="histogram", fill=I("#5ab4ac"),
      col=I("white"), binwidth=1) + ggtitle("Instagram (n = 128)") + xlab("Year") +
  ylab("No. Occurances")

# For iNaturalist alone:
#nrow(moth[moth$source == "iNaturalist",]) # 1442
#qplot(moth$year[moth$source=="iNaturalist"], geom="histogram", fill=I("#5ab4ac"),
#      col=I("white"), binwidth=1) + ggtitle("iNaturalist (n = 512)") + xlab("Year") +
#  ylab("No. Occurances")

#### Showing how distribution has changed over time ####
## Build map with only nations of interest
data("wrld_simpl") #Get border data from rgeos
myn <- c("United Kingdom", "France","Italy", "Germany", "Austria", "Germany", "Denmark",
         "Ireland", "Switzerland", "Belgium", "Netherlands", "Luxembourg", "Czech Republic")
my_map <- wrld_simpl[wrld_simpl$NAME %in% myn,]
plot(my_map, axes=T)

## And for ggspatial

world <- ne_countries(scale = "medium", returnclass = "sf")
myn2 <- c("United Kingdom", "France","Italy", "Germany", "Austria", "Germany", "Denmark",
          "Ireland", "Switzerland", "Belgium", "Netherlands", "Luxembourg", "Czech Rep.")
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
#plot(wrld_simpl, ylim=c(35, 60), xlim=c(-10, 35), axes=T,
#     main="iNaturalist")
#points(moth$Long[moth$source=="iNaturalist"],
#       moth$Lat[moth$source=="iNaturalist"], col="dark green")
#points(moth$Long[moth$source=="GBIF"],
#       moth$Lat[moth$source=="GBIF"], col="red")
