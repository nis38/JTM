###### Boxplots

## Includes extra sensitivity analyses

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
  geom_signif(comparisons = list(c("Flickr", "GBIF")), map_signif_level = T) +
  geom_signif(comparisons = list(c("Instagram", "Flickr")), map_signif_level = T) +
  geom_signif(comparisons = list(c("GBIF", "Instagram")), map_signif_level = T)

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
  geom_signif(comparisons = list(c("Flickr", "GBIF")), map_signif_level = T) +
  geom_signif(comparisons = list(c("Flickr", "Instagram")), map_signif_level = T) +
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