#######################
######### HSMs #######
######################

#### 00s (Baseline Model) ####
# climate input
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
  
  res <- list(hsm, bc.eval@auc, predict.presence)
  return(res)
  
}

occur00s <- as.data.frame(cbind(moth$Long[moth$year > 2000 & moth$year < 2009 & moth$source=="GBIF"],
                                moth$Lat[moth$year > 2000 & moth$year < 2009 & moth$source=="GBIF"]))


hsmodplot <- function(x, y, a, z) {
  colnames(y) <- c("Long", "Lat")
  colnames(x) <- c("Long", "Lat")
  ggplot(data = my_mapg) +
    geom_sf(color = "black", fill = "white") +
    ylim(c(35, 60)) +
    xlim(c(-10, 20)) +
    layer_spatial(z > bc.threshold) +
    scale_fill_binned(low = "white", high = "#999933", na.value = "transparent") +
    xlab("Longitude") +
    ylab("Latitude") +
    geom_sf(color = "black", fill = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(color=c(1,0,1,0)),
          axis.text.y = element_text(color=c(1,0,1,0))) +
    geom_point(data = x, mapping = aes(x = Long, y = Lat), color = "black", size = 1) +
    geom_point(data = y, mapping = aes(x = Long, y = Lat), color = "turquoise", size = 1) +
    geom_point(data = a, mapping = aes(x = Long, y = Lat), color = "red", size = 1) +
    # change the number of point lines depending on data input - needed for 2010, 2011, 2013 data because no instagram records
    theme(legend.position = "none")
}
## 10s models ####

prec2010 <- raster("Climate/2010/prec2010.tif")
covprec2010 <- raster("Climate/2010/covprec2010.tif")
tmax2010 <- raster("Climate/2010/tmax2010.tif")
covtmax2010 <- raster("Climate/2010/covtmax2010.tif")
nl2010 <- raster("Climate/2010/nl2012.tif")

clim2010 <- stack(prec2010, covprec2010, tmax2010, covtmax2010)
clim2010 <- mask(clim2010, my_map)
occur2010 <- as.data.frame(cbind(moth$Long[moth$year == 2010 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2010 & moth$source=="GBIF"]))

hsm2010 <- hsmmod(clim2010, occur00s)
hsm2010 <- hsm2010[[1]]
plot(hsm2010)

hsmauc2010 <- replicate(n = 5, hsmmod(clim2010, occur00s), simplify = F)
hsmauc2010m <- mean(c(hsmauc2010[[1]][[2]], hsmauc2010[[2]][[2]],
                      hsmauc2010[[3]][[2]], hsmauc2010[[4]][[2]],
                      hsmauc2010[[5]][[2]]))
hsmauc2010m

ioccur2010 <- as.data.frame(cbind(moth$Long[moth$year == 2010 & moth$source!="GBIF"],
                                  moth$Lat[moth$year == 2010 & moth$source!="GBIF"]))

hsmodplot(occur2010, ioccur2010, hsmauc2010[[1]][[3]])

prec2011 <- raster("Climate/2011/prec2011.tif")
covprec2011 <- raster("Climate/2011/covprec2011.tif")
tmax2011 <- raster("Climate/2011/tmax2011.tif")
covtmax2011 <- raster("Climate/2011/covtmax2011.tif")
nl2011 <- raster("Climate/2011/nl2012.tif")

clim2011 <- stack(prec2011, covprec2011, tmax2011, covtmax2011)
clim2011 <- mask(clim2011, my_map)
occur2011 <- as.data.frame(cbind(moth$Long[moth$year == 2011 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2011 & moth$source=="GBIF"]))

hsmauc2011 <- replicate(n = 5, hsmmod(clim2011, occur00s), simplify = F)
hsmauc2011m <- mean(c(hsmauc2011[[1]][[2]], hsmauc2011[[2]][[2]], hsmauc2011[[3]][[2]],
                      hsmauc2011[[4]][[2]], hsmauc2011[[5]][[2]]))
hsmauc2011m

ioccur2011 <- as.data.frame(cbind(moth$Long[moth$year == 2011 & moth$source=="Instagram"],
                                  moth$Lat[moth$year == 2011 & moth$source=="Instagram"]))
foccur2011 <- as.data.frame(cbind(moth$Long[moth$year == 2011 & moth$source=="Flickr"],
                                  moth$Lat[moth$year == 2011 & moth$source=="Flickr"]))

hsmodplot(occur2011, ioccur2011, foccur2011, hsmauc2011[[1]][[3]])


prec2012 <- raster("Climate/2012/prec2012.tif")
covprec2012 <- raster("Climate/2012/covprec2012.tif")
tmax2012 <- raster("Climate/2012/tmax2012.tif")
covtmax2012 <- raster("Climate/2012/covtmax2012.tif")
nl2012 <- raster("Climate/2012/nl2012.tif")

clim2012 <- stack(prec2012, covprec2012, tmax2012, covtmax2012)
clim2012 <- mask(clim2012, my_map)
occur2012 <- as.data.frame(cbind(moth$Long[moth$year == 2012 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2012 & moth$source=="GBIF"]))

hsmauc2012 <- replicate(n = 5, hsmmod(clim2012, occur00s), simplify = F)
hsmauc2012m <- mean(c(hsmauc2012[[1]][[2]], hsmauc2012[[2]][[2]],
                      hsmauc2012[[3]][[2]], hsmauc2012[[4]][[2]], hsmauc2012[[5]][[2]]))
hsmauc2012m

ioccur2012 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2012 & moth$source=="Instagram"],
                                  Lat = moth$Lat[moth$year == 2012 & moth$source=="Instagram"]))
foccur2012 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2012 & moth$source=="Flickr"],
                                  Lat = moth$Lat[moth$year == 2012 & moth$source=="Flickr"]))

hsmodplot(occur2012, ioccur2012, foccur2012, hsmauc2012[[1]][[3]])


prec2013 <- raster("Climate/2013/prec2013.tif")
covprec2013 <- raster("Climate/2013/covprec2013.tif")
tmax2013 <- raster("Climate/2013/tmax2013.tif")
covtmax2013 <- raster("Climate/2013/covtmax2013.tif")
nl2013 <- raster("Climate/2013/nl2013.tif")

clim2013 <- stack(prec2013, covprec2013, tmax2013, covtmax2013)
clim2013 <- mask(clim2013, my_map)
occur2013 <- as.data.frame(cbind(moth$Long[moth$year == 2013 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2013 & moth$source=="GBIF"]))

hsmauc2013 <- replicate(n = 5, hsmmod(clim2013, occur00s), simplify = F)
hsmauc2013m <- mean(c(hsmauc2013[[1]][[2]], hsmauc2013[[2]][[2]], hsmauc2013[[3]][[2]],
                      hsmauc2013[[4]][[2]], hsmauc2013[[5]][[2]]))
hsmauc2013m

ioccur2013 <- as.data.frame(cbind(moth$Long[moth$year == 2013 & moth$source!="GBIF"],
                                  moth$Lat[moth$year == 2013 & moth$source!="GBIF"]))

hsmodplot(occur2013, ioccur2013, hsmauc2013[[1]][[3]])

prec2014 <- raster("Climate/2014/prec2014.tif")
covprec2014 <- raster("Climate/2014/covprec2014.tif")
tmax2014 <- raster("Climate/2014/tmax2014.tif")
covtmax2014 <- raster("Climate/2014/covtmax2014.tif")
nl2014 <- raster("Climate/2014/nl2014.tif")

clim2014 <- stack(prec2014, covprec2014, tmax2014, covtmax2014)
clim2014 <- mask(clim2014, my_map)
occur2014 <- as.data.frame(cbind(moth$Long[moth$year == 2014 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2014 & moth$source=="GBIF"]))

hsmauc2014 <- replicate(n = 5, hsmmod(clim2014, occur00s), simplify = F)
hsmauc2014m <- mean(c(hsmauc2014[[1]][[2]], hsmauc2014[[2]][[2]], hsmauc2014[[3]][[2]],
                      hsmauc2014[[4]][[2]], hsmauc2014[[5]][[2]]))
hsmauc2014m

ioccur2014 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2014 & moth$source=="Instagram"],
                                  Lat = moth$Lat[moth$year == 2014 & moth$source=="Instagram"]))
foccur2014 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2014 & moth$source=="Flickr"],
                                  Lat = moth$Lat[moth$year == 2014 & moth$source=="Flickr"]))

hsmodplot(occur2014, ioccur2014, foccur2014, hsmauc2014[[1]][[3]])


prec2015 <- raster("Climate/2015/prec2015.tif")
covprec2015 <- raster("Climate/2015/covprec2015.tif")
tmax2015 <- raster("Climate/2015/tmax2015.tif")
covtmax2015 <- raster("Climate/2015/covtmax2015.tif")
nl2015 <- raster("Climate/2015/nl2015.tif")

clim2015 <- stack(prec2015, covprec2015, tmax2015, covtmax2015)
clim2015 <- mask(clim2015, my_map)
occur2015 <- as.data.frame(cbind(moth$Long[moth$year == 2015 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2015 & moth$source=="GBIF"]))

hsmauc2015 <- replicate(n = 5, hsmmod(clim2015, occur00s), simplify = F)
hsmauc2015m <- mean(c(hsmauc2015[[1]][[2]], hsmauc2015[[2]][[2]], hsmauc2015[[3]][[2]],
                      hsmauc2015[[4]][[2]], hsmauc2015[[5]][[2]]))
hsmauc2015m

ioccur2015 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2015 & moth$source=="Instagram"],
                                  Lat = moth$Lat[moth$year == 2015 & moth$source=="Instagram"]))
foccur2015 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2015 & moth$source=="Flickr"],
                                  Lat = moth$Lat[moth$year == 2015 & moth$source=="Flickr"]))


hsmodplot(occur2015, ioccur2015, foccur2015, hsmauc2015[[1]][[3]])


prec2016 <- raster("Climate/2016/prec2016.tif")
covprec2016 <- raster("Climate/2016/covprec2016.tif")
tmax2016 <- raster("Climate/2016/tmax2016.tif")
covtmax2016 <- raster("Climate/2016/covtmax2016.tif")
nl2016 <- raster("Climate/2016/nl2016.tif")

clim2016 <- stack(prec2016, covprec2016, tmax2016, covtmax2016)
clim2016 <- mask(clim2016, my_map)
occur2016 <- as.data.frame(cbind(moth$Long[moth$year == 2016 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2016 & moth$source=="GBIF"]))

hsmauc2016 <- replicate(n = 5, hsmmod(clim2016, occur00s), simplify = F)
hsmauc2016m <- mean(c(hsmauc2016[[1]][[2]], hsmauc2016[[2]][[2]], hsmauc2016[[3]][[2]],
                      hsmauc2016[[4]][[2]], hsmauc2016[[5]][[2]]))
hsmauc2016m

ioccur2016 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2016 & moth$source=="Instagram"],
                                  Lat = moth$Lat[moth$year == 2016 & moth$source=="Instagram"]))
foccur2016 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2016 & moth$source=="Flickr"],
                                  Lat = moth$Lat[moth$year == 2016 & moth$source=="Flickr"]))


hsmodplot(occur2016, ioccur2016, foccur2016, hsmauc2016[[1]][[3]])
hsm16 <- hsm


prec2017 <- raster("Climate/2017/prec2017.tif")
covprec2017 <- raster("Climate/2017/covprec2017.tif")
tmax2017 <- raster("Climate/2017/tmax2017.tif")
covtmax2017 <- raster("Climate/2017/covtmax2017.tif")
nl2017 <- raster("Climate/2017/nl2017.tif")

clim2017 <- stack(prec2017, covprec2017, tmax2017, covtmax2017)
clim2017 <- mask(clim2017, my_map)
occur2017 <- as.data.frame(cbind(moth$Long[moth$year == 2017 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2017 & moth$source=="GBIF"]))

hsmauc2017 <- replicate(n = 5, hsmmod(clim2017, occur00s), simplify = F)
hsmauc2017m <- mean(c(hsmauc2017[[1]][[2]], hsmauc2017[[2]][[2]], hsmauc2017[[3]][[2]],
                      hsmauc2017[[4]][[2]], hsmauc2017[[5]][[2]]))
hsmauc2017m

ioccur2017 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2017 & moth$source=="Instagram"],
                                  Lat = moth$Lat[moth$year == 2017 & moth$source=="Instagram"]))
foccur2017 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2017 & moth$source=="Flickr"],
                                  Lat = moth$Lat[moth$year == 2017 & moth$source=="Flickr"]))


hsmodplot(occur2017, ioccur2017, foccur2017, hsmauc2017[[1]][[3]])
hsm17 <- hsm


prec2018 <- raster("Climate/2018/prec2018.tif")
covprec2018 <- raster("Climate/2018/covprec2018.tif")
tmax2018 <- raster("Climate/2018/tmax2018.tif")
covtmax2018 <- raster("Climate/2018/covtmax2018.tif")
nl2018 <- raster("Climate/2018/nl2018.tif")

clim2018 <- stack(prec2018, covprec2018, tmax2018, covtmax2018)
clim2018 <- mask(clim2018, my_map)
occur2018 <- as.data.frame(cbind(moth$Long[moth$year == 2018 & moth$source=="GBIF"],
                                 moth$Lat[moth$year == 2018 & moth$source=="GBIF"]))

hsmauc2018 <- replicate(n = 5, hsmmod(clim2018, occur00s), simplify = F)
hsmauc2018m <- mean(c(hsmauc2018[[1]][[2]], hsmauc2018[[2]][[2]], hsmauc2018[[3]][[2]],
                      hsmauc2018[[4]][[2]], hsmauc2018[[5]][[2]]))
hsmauc2018m

ioccur2018 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2018 & moth$source=="Instagram"],
                                  Lat = moth$Lat[moth$year == 2018 & moth$source=="Instagram"]))
foccur2018 <- as.data.frame(cbind(Long = moth$Long[moth$year == 2018 & moth$source=="Flickr"],
                                  Lat = moth$Lat[moth$year == 2018 & moth$source=="Flickr"]))


hsmodplot(occur2018, ioccur2018, foccur2018, hsmauc2018[[1]][[3]])
hsm18 <- hsm