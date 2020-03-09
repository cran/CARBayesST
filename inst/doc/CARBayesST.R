### R code from vignette source 'CARBayesST.Rnw'

###################################################
### code chunk number 1: CARBayesST.Rnw:92-93
###################################################
options(prompt = "R>  ")


###################################################
### code chunk number 2: CARBayesST.Rnw:429-433
###################################################
library("CARBayesdata")
library("sp")
data("GGHB.IG")
data("pollutionhealthdata")


###################################################
### code chunk number 3: CARBayesST.Rnw:439-440
###################################################
head(pollutionhealthdata)


###################################################
### code chunk number 4: CARBayesST.Rnw:445-449
###################################################
library(dplyr)
pollutionhealthdata <- pollutionhealthdata %>% mutate( 
          SMR = pollutionhealthdata$observed / pollutionhealthdata$expected, 
          logSMR = log(pollutionhealthdata$observed / pollutionhealthdata$expected))


###################################################
### code chunk number 5: CARBayesST.Rnw:455-457
###################################################
library(GGally)
ggpairs(pollutionhealthdata, columns=c(9, 5:7))


###################################################
### code chunk number 6: CARBayesST.Rnw:472-475
###################################################
group_IG <- group_by(pollutionhealthdata, IG)
SMR.av <- summarise(group_IG, SMR.mean = mean(SMR))
GGHB.IG@data$SMR <- SMR.av$SMR.mean


###################################################
### code chunk number 7: CARBayesST.Rnw:482-484
###################################################
library(rgdal)
GGHB.IG <- spTransform(GGHB.IG, CRS("+proj=longlat +datum=WGS84 +no_defs"))


###################################################
### code chunk number 8: CARBayesST.Rnw:489-499
###################################################
library(leaflet)
colours <- colorNumeric(palette = "YlOrRd", domain = GGHB.IG@data$SMR)
map1 <- leaflet(data=GGHB.IG) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(SMR), weight=1, color="",
                fillOpacity = 0.7) %>%
    addLegend(pal = colours, values = GGHB.IG@data$SMR, opacity = 1, 
                title="SMR") %>%
    addScaleBar(position="bottomleft")
map1


###################################################
### code chunk number 9: CARBayesST.Rnw:514-518
###################################################
library("spdep")
W.nb <- poly2nb(GGHB.IG, row.names = SMR.av$IG)
W.list <- nb2listw(W.nb, style = "B")
W <- nb2mat(W.nb, style = "B")


###################################################
### code chunk number 10: CARBayesST.Rnw:527-533
###################################################
formula <- observed ~ offset(log(expected)) + jsa + price + pm10
model1 <- glm(formula = formula, family = "quasipoisson", 
    data = pollutionhealthdata)
resid.glm <- residuals(model1)
summary(model1)$coefficients
summary(model1)$dispersion


###################################################
### code chunk number 11: CARBayesST.Rnw:538-539
###################################################
moran.mc(x = resid.glm[1:271], listw = W.list, nsim = 10000)


###################################################
### code chunk number 12: CARBayesST.Rnw:692-697
###################################################
library("CARBayesdata")
library("sp")
data("GGHB.IG")
data("salesdata")
head(salesdata)


###################################################
### code chunk number 13: CARBayesST.Rnw:703-710
###################################################
salesdata <- salesdata %>% mutate(salesprop = salesdata$sales / salesdata$stock)
library(ggplot2)
ggplot(salesdata, aes(x = factor(year), y = salesprop)) +
    geom_boxplot(fill="red", alpha=0.7) + 
    scale_x_discrete(name = "Year") +
    scale_y_continuous(name = "Sales proportion") + 
    theme(text=element_text(size=16), plot.title=element_text(size=18, face="bold")) 


###################################################
### code chunk number 14: CARBayesST.Rnw:723-727
###################################################
library(dplyr)
group_IG <- group_by(salesdata, IG)
salesprop <- summarise(group_IG, salesproprtion.mean = mean(salesprop))
GGHB.IG@data$sales <- salesprop$salesproprtion.mean


###################################################
### code chunk number 15: CARBayesST.Rnw:732-745
###################################################
library(rgdal)
GGHB.IG <- spTransform(GGHB.IG, CRS("+proj=longlat +datum=WGS84 +no_defs"))

library(leaflet)
colours <- colorNumeric(palette = "YlOrRd", domain = GGHB.IG@data$sales)
map1 <- leaflet(data=GGHB.IG) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(sales), color="", weight=1, 
                fillOpacity = 0.7) %>%
    addLegend(pal = colours, values = GGHB.IG@data$sales, opacity = 1, 
                title="Sales") %>%
    addScaleBar(position="bottomleft")
map1


###################################################
### code chunk number 16: CARBayesST.Rnw:763-766
###################################################
library("spdep")
W.nb <- poly2nb(GGHB.IG, row.names = salesprop$salesproprtion.mean)
W <- nb2mat(W.nb, style = "B")


