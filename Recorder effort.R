#####################
## Recorder effort ##
#####################

## requires data from EBBA2.

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
#bbi <- read.csv("bbREiNat.csv", header = T, stringsAsFactors = F)
#colnames(bbi) [c(1, 3, 22, 23)] <- c("ID", "year", "Lat", "Long")
#keep5 <- c("ID", "year", "Lat", "Long")
#bbi <- bbi[keep5]
#bbi$year <- substr(bbi$year, 0, 4)
#bbi$Source <- "iNaturalist"
#bbi <- bbi[bbi$year>1999 & bbi$year<2019,]
#bbi <- bbi[complete.cases(bbi$Lat),]
#bbi <- bbi[complete.cases(bbi$Long),]

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


#bbi2016 <- bbrast(bbi[bbi$year == 2016,])
#bbi2016[is.na(bbi2016)] <- 0

#bb@data$RE.i16 <- raster::extract(bbi2016, 
#                                  bb,
#                                  fun = sum)
#bb@data$RE.i16 <- bb@data$RE.i16/bb@data$abundance_2

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
keep9 <- c("HS", "RE.G16", "RE.F16", "RE.I16")
REdf16 <- REdf16[keep9]
REdf16$RE.G16 <- ifelse(REdf16$RE.G16==Inf, 0, REdf16$RE.G16)
#REdf16$RE.i16 <- ifelse(REdf16$RE.i16==Inf, 0, REdf16$RE.i16)
REdf16$RE.F16 <- ifelse(REdf16$RE.F16==Inf, 0, REdf16$RE.F16)
REdf16$RE.I16 <- ifelse(REdf16$RE.I16==Inf, 0, REdf16$RE.I16)

RE2016 <- stack(hsm16, bbG2016, bbF2016, bbI2016)

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
#rp16$RE.i16 <- ifelse(rp16$RE.i16==Inf, 0, rp16$RE.i16)
rp16$RE.F16 <- ifelse(rp16$RE.F16==Inf, 0, rp16$RE.F16)
rp16$RE.I16 <- ifelse(rp16$RE.I16==Inf, 0, rp16$RE.I16)

sq16 <- rbind(REdf16, rp16)
string16 <- c(rep(1, nrow(REdf16)), rep(0, nrow(rp16)))

REdf16 <- data.frame(cbind(pa=string16, sq16))
colnames(REdf16)

occurb2016 <- as.data.frame(cbind(moth$source[moth$year==2016]))
occurb2016 <- rbind(occurb2016, occurb2016)
REdf16 <- cbind(REdf16, occurb2016)
colnames(REdf16) [6] <- "source"

##

bbG2017 <- bbrast(bbG[bbG$year == 2017,])
bbG2017[is.na(bbG2017)] <- 0

bb@data$RE.G17 <- raster::extract(bbG2017, 
                                  bb,
                                  fun = sum)
bb@data$RE.G17 <- bb@data$RE.G17/bb@data$abundance_2


#bbi2017 <- bbrast(bbi[bbi$year == 2017,])
#bbi2017[is.na(bbi2017)] <- 0

#bb@data$RE.i17 <- raster::extract(bbi2017, 
#                                  bb,
#                                  fun = sum)
#bb@data$RE.i17 <- bb@data$RE.i17/bb@data$abundance_2

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
keep9 <- c("HS", "RE.G17", "RE.F17", "RE.I17")
REdf17 <- REdf17[keep9]
REdf17$RE.G17 <- ifelse(REdf17$RE.G17==Inf, 0, REdf17$RE.G17)
#REdf17$RE.i17 <- ifelse(REdf17$RE.i17==Inf, 0, REdf17$RE.i17)
REdf17$RE.F17 <- ifelse(REdf17$RE.F17==Inf, 0, REdf17$RE.F17)
REdf17$RE.I17 <- ifelse(REdf17$RE.I17==Inf, 0, REdf17$RE.I17)

RE2017 <- stack(hsm17, bbG2017, bbF2017, bbI2017)
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
#rp17$RE.i17 <- ifelse(rp17$RE.i17==Inf, 0, rp17$RE.i17)
rp17$RE.F17 <- ifelse(rp17$RE.F17==Inf, 0, rp17$RE.F17)
rp17$RE.I17 <- ifelse(rp17$RE.I17==Inf, 0, rp17$RE.I17)

sq17 <- rbind(REdf17, rp17) 
string17 <- c(rep(1, nrow(REdf17)), rep(0, nrow(rp17)))

REdf17 <- data.frame(cbind(pa=string17, sq17))
colnames(REdf17)

occurb2017 <- as.data.frame(cbind(moth$source[moth$year==2017]))
occurb2017 <- rbind(occurb2017, occurb2017)
REdf17 <- cbind(REdf17, occurb2017)
colnames(REdf17) [6] <- "source"

##

bbG2018 <- bbrast(bbG[bbG$year == 2018,])
bbG2018[is.na(bbG2018)] <- 0

bb@data$RE.G18 <- raster::extract(bbG2018, 
                                  bb,
                                  fun = sum)
bb@data$RE.G18 <- bb@data$RE.G18/bb@data$abundance_2


#bbi2018 <- bbrast(bbi[bbi$year == 2018,])
#bbi2018[is.na(bbi2018)] <- 0

#bb@data$RE.i18 <- raster::extract(bbi2018, 
#                                  bb,
#                                  fun = sum)
#bb@data$RE.i18 <- bb@data$RE.i18/bb@data$abundance_2

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
keep9 <- c("HS", "RE.G18", "RE.F18", "RE.I18")
REdf18 <- REdf18[keep9]
REdf18$RE.G18 <- ifelse(REdf18$RE.G18==Inf, 0, REdf18$RE.G18)
#REdf18$RE.i18 <- ifelse(REdf18$RE.i18==Inf, 0, REdf18$RE.i18)
REdf18$RE.F18 <- ifelse(REdf18$RE.F18==Inf, 0, REdf18$RE.F18)
REdf18$RE.I18 <- ifelse(REdf18$RE.I18==Inf, 0, REdf18$RE.I18)

RE2018 <- stack(hsm18, bbG2018, bbF2018, bbI2018)

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
#rp18$RE.i18 <- ifelse(rp18$RE.i18==Inf, 0, rp18$RE.i18)
rp18$RE.F18 <- ifelse(rp18$RE.F18==Inf, 0, rp18$RE.F18)
rp18$RE.I18 <- ifelse(rp18$RE.I18==Inf, 0, rp18$RE.I18)

sq18 <- rbind(REdf18, rp18) 
string18 <- c(rep(1, nrow(REdf18)), rep(0, nrow(rp18)))

REdf18 <- data.frame(cbind(pa=string18, sq18))
colnames(REdf18)

occurb2018 <- as.data.frame(cbind(moth$source[moth$year==2018]))
occurb2018 <- rbind(occurb2018, occurb2018)
REdf18 <- cbind(REdf18, occurb2018)
colnames(REdf18) [6] <- "source"

##

colnames(REdf16) [c(3, 4, 5, 6)] <- c("RE.G", "RE.F", "RE.I", "source")
colnames(REdf17) [c(3, 4, 5)] <- c("RE.G", "RE.F", "RE.I")
colnames(REdf18) [c(3, 4, 5, 6)] <- c("RE.G", "RE.F", "RE.I",
                                      "source")

REdf <- rbind(REdf16, REdf17, REdf18)
REdf$RE.G <- ifelse(is.na(REdf$RE.G), 0, REdf$RE.G)
REdf$RE.F <- ifelse(is.na(REdf$RE.F), 0, REdf$RE.F)
REdf$RE.I <- ifelse(is.na(REdf$RE.I), 0, REdf$RE.I)

REdfG <- REdf[REdf$source=="GBIF",]
REdfF <- REdf[REdf$source=="Flickr",]
REdfI <- REdf[REdf$source=="Instagram",]

REdfG$mRE <- mean(as.numeric(REdfG$RE.G))
REdfG$sdRE <- sd(as.numeric(REdfG$RE.G))
REdfG$sRE <- (as.numeric(REdfG$RE.G)-as.numeric(REdfG$mRE))/as.numeric(REdfG$sdRE)

#REdfi$mRE <- mean(REdfi$RE.i)
#REdfi$sdRE <- sd(REdfi$RE.i)
#REdfi$sRE <- (REdfi$RE.i-REdfi$mRE)/REdfi$sdRE

REdfF$mRE <- mean(as.numeric(REdfF$RE.G))
REdfF$sdRE <- sd(as.numeric(REdfF$RE.G))
REdfF$sRE <- (as.numeric(REdfF$RE.G)-as.numeric(REdfF$mRE))/as.numeric(REdfF$sdRE)

REdfI$mRE <- mean(as.numeric(REdfI$RE.G))
REdfI$sdRE <- sd(as.numeric(REdfI$RE.G))
REdfI$sRE <- (as.numeric(REdfI$RE.G)-as.numeric(REdfI$mRE))/as.numeric(REdfI$sdRE)

REdfG <- REdfG[complete.cases(REdfG$HS),]
#REdfi <- REdfi[complete.cases(REdfi$HS),]
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
#glmmi <- glmmTMB(pa ~ HS + sRE,
#                 family=binomial(link = "logit"), data = REdfi)
#dredge(glmmi)#
#
#ggHSi <- plot(ggpredict(glmmi, c("HS [all]")), add = T,
#              dot.size = 0.3, show.title = F) +
#  labs(x = "Habitat Suitability", y = "Presence/Absence", col = "Year") +
#  scale_y_continuous(breaks = c(0, 1)) + theme_bw() + ggtitle("iNaturalist")
#
#ggREi <- plot(ggpredict(glmmi, c("sRE [all]")), add = T,
#              dot.size = 0.7, show.title = F) +
#  labs(x = "Recorder Effort", y = "Presence/Absence", col = "Year") +
#  scale_y_continuous(breaks = c(0, 1)) + xlim(c(-0.5, 3)) + theme_bw()
#
#plot_grid(ggHSi, ggREi, nrow = 2, labels = c("(a)", "(b)"),
#          label_size = 10)
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