### R code from vignette source 'CARBayesST.Rnw'

###################################################
### code chunk number 1: CARBayesST.Rnw:92-93
###################################################
options(prompt = "R>  ")


###################################################
### code chunk number 2: CARBayesST.Rnw:537-541
###################################################
library("CARBayesdata")
library("sf")
data("GGHB.IZ")
data("pollutionhealthdata")


###################################################
### code chunk number 3: CARBayesST.Rnw:547-551
###################################################
class(GGHB.IZ)
head(GGHB.IZ)
class(pollutionhealthdata)
head(pollutionhealthdata)


###################################################
### code chunk number 4: CARBayesST.Rnw:556-561
###################################################
library(dplyr)
pollutionhealthdata <- pollutionhealthdata %>% mutate( 
          SMR = pollutionhealthdata$observed / pollutionhealthdata$expected, 
          logSMR = log(pollutionhealthdata$observed / pollutionhealthdata$expected))
head(pollutionhealthdata)


###################################################
### code chunk number 5: CARBayesST.Rnw:567-569
###################################################
library(GGally)
ggpairs(pollutionhealthdata, columns=c(9, 5:7))


###################################################
### code chunk number 6: CARBayesST.Rnw:584-588
###################################################
group_IZ <- group_by(pollutionhealthdata, IZ)
SMR.av <- summarise(group_IZ, SMR.mean = mean(SMR))
GGHB.IZ$SMR <- SMR.av$SMR.mean
head(GGHB.IZ)


###################################################
### code chunk number 7: CARBayesST.Rnw:595-596
###################################################
GGHB.IZ <- st_transform(x=GGHB.IZ, crs='+proj=longlat +datum=WGS84 +no_defs')


###################################################
### code chunk number 8: CARBayesST.Rnw:601-611
###################################################
library(leaflet)
colours <- colorNumeric(palette = "YlOrRd", domain = GGHB.IZ$SMR)
leaflet(data=GGHB.IZ) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(GGHB.IZ$SMR), 
                color="grey", weight=1, 
                fillOpacity = 0.7) %>%
    addLegend(pal = colours, values = GGHB.IZ$SMR, 
              opacity = 1, title="SMR") %>%
    addScaleBar(position="bottomleft")


###################################################
### code chunk number 9: CARBayesST.Rnw:626-630
###################################################
library("spdep")
W.nb <- poly2nb(GGHB.IZ, row.names = GGHB.IZ$IZ)
W.list <- nb2listw(W.nb, style = "B")
W <- nb2mat(W.nb, style = "B")


###################################################
### code chunk number 10: CARBayesST.Rnw:639-645
###################################################
formula <- observed ~ offset(log(expected)) + jsa + price + pm10
model1 <- glm(formula = formula, family = "quasipoisson", 
    data = pollutionhealthdata)
resid.glm <- residuals(model1)
summary(model1)$coefficients
summary(model1)$dispersion


###################################################
### code chunk number 11: CARBayesST.Rnw:650-651
###################################################
moran.mc(x = resid.glm[1:271], listw = W.list, nsim = 10000)


###################################################
### code chunk number 12: CARBayesST.Rnw:783-788
###################################################
library("CARBayesdata")
library("sf")
data("GGHB.IZ")
data("salesdata")
head(salesdata)


###################################################
### code chunk number 13: CARBayesST.Rnw:794-801
###################################################
salesdata <- salesdata %>% mutate(salesprop = salesdata$sales / salesdata$stock)
library(ggplot2)
ggplot(salesdata, aes(x = factor(year), y = salesprop)) +
    geom_boxplot(fill="red", alpha=0.7) + 
    scale_x_discrete(name = "Year") +
    scale_y_continuous(name = "Sales proportion") + 
    theme(text=element_text(size=16), plot.title=element_text(size=18, face="bold")) 


###################################################
### code chunk number 14: CARBayesST.Rnw:814-819
###################################################
library(dplyr)
group_IZ <- group_by(salesdata, IZ)
salesprop <- summarise(group_IZ, salesproprtion.mean = mean(salesprop))
GGHB.IZ$sales <- salesprop$salesproprtion.mean
head(GGHB.IZ)


###################################################
### code chunk number 15: CARBayesST.Rnw:824-836
###################################################
GGHB.IZ <- st_transform(x=GGHB.IZ, crs='+proj=longlat +datum=WGS84 +no_defs')

library(leaflet)
colours <- colorNumeric(palette = "YlOrRd", domain = GGHB.IZ$sales)
map1 <- leaflet(data=GGHB.IZ) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(sales), color="grey", weight=1, 
                fillOpacity = 0.7) %>%
    addLegend(pal = colours, values = GGHB.IZ$sales, opacity = 1, 
              title="Sales") %>%
    addScaleBar(position="bottomleft")
map1


###################################################
### code chunk number 16: CARBayesST.Rnw:854-857
###################################################
library("spdep")
W.nb <- poly2nb(GGHB.IZ, row.names = GGHB.IZ$IZ)
W <- nb2mat(W.nb, style = "B")


