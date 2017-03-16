### R code from vignette source 'CARBayesST.Rnw'

###################################################
### code chunk number 1: CARBayesST.Rnw:95-96
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 2: CARBayesST.Rnw:378-385
###################################################
n.space <- 20
N <- 20
x.easting <- 1:n.space
x.northing <- 1:n.space
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)
N.all <- N * K


###################################################
### code chunk number 3: CARBayesST.Rnw:390-393
###################################################
distance <- as.matrix(dist(Grid))
W <- array(0, c(K,K))
W[distance==1] <- 1


###################################################
### code chunk number 4: CARBayesST.Rnw:405-408
###################################################
distance <- as.matrix(dist(1:N))
D <-array(0, c(N,N))
D[distance==1] <-1


###################################################
### code chunk number 5: CARBayesST.Rnw:413-414
###################################################
Q.W <- 0.8 * (diag(apply(W, 2, sum)) - W) + 0.2 * diag(rep(1,K))


###################################################
### code chunk number 6: CARBayesST.Rnw:419-423
###################################################
Q.W.inv <- solve(Q.W)
library("MASS")
phi <- mvrnorm(n = 1, mu = rep(0, K), Sigma = (0.01 * Q.W.inv))
phi.long <- rep(phi, N)


###################################################
### code chunk number 7: CARBayesST.Rnw:429-433
###################################################
Q.D <- 0.8 * (diag(apply(D, 2, sum)) - D) + 0.2 * diag(rep(1, N))
Q.D.inv <- solve(Q.D)
delta <- mvrnorm(n = 1, mu = rep(0, N), Sigma = (0.01 * Q.D.inv))
delta.long <- kronecker(delta, rep(1, K))


###################################################
### code chunk number 8: CARBayesST.Rnw:438-440
###################################################
x <- rnorm(n = N.all, mean = 0, sd = 1)
gamma <- rnorm(n = N.all, mean = 0, sd = sqrt(0.01))


###################################################
### code chunk number 9: CARBayesST.Rnw:445-451
###################################################
beta1 <- 0
beta2 <- 0.1
n <- rep(50, N.all)
LP <- beta1 + beta2 * x + phi.long + delta.long + gamma
theta.true <- exp(LP) / (1 + exp(LP))
Y <- rbinom(n = N.all, size = n, prob = theta.true)


###################################################
### code chunk number 10: CARBayesST.Rnw:587-591
###################################################
library("CARBayesdata")
library("sp")
data("GGHB.IG")
data("pollutionhealthdata")


###################################################
### code chunk number 11: CARBayesST.Rnw:597-598
###################################################
head(pollutionhealthdata)


###################################################
### code chunk number 12: CARBayesST.Rnw:604-609
###################################################
pollutionhealthdata$SMR <- pollutionhealthdata$observed / pollutionhealthdata$expected
pollutionhealthdata$logSMR <- log(pollutionhealthdata$SMR)
par(pty="s", cex.axis=1.5, cex.lab=1.5)
pairs(pollutionhealthdata[ ,c(9, 5:7)], pch=19, cex=0.5, lower.panel=NULL, panel=panel.smooth,
      labels=c("ln(SMR)", "PM10", "JSA", "Price (*100,000)"))


###################################################
### code chunk number 13: CARBayesST.Rnw:633-637
###################################################
library("dplyr")
SMR.av <- summarise(group_by(pollutionhealthdata,IG), SMR.mean = 
    mean(SMR))
GGHB.IG@data$SMR <- SMR.av$SMR.mean


###################################################
### code chunk number 14: CARBayesST.Rnw:643-656
###################################################
l1 = list("SpatialPolygonsRescale", layout.north.arrow(), 
    offset = c(220000,647000), scale = 4000)
l2 = list("SpatialPolygonsRescale", layout.scale.bar(), offset = 
    c(225000, 647000), scale = 10000, fill = c("transparent","black"))
l3 = list("sp.text", c(225000,649000), "0")
l4 = list("sp.text", c(230000,649000), "5000 m")
breakpoints <- seq(min(SMR.av$SMR.mean)-0.1, max(SMR.av$SMR.mean)+0.1, 
    length.out = 11)
spplot(GGHB.IG, "SMR", sp.layout = list(l1, l2, l3, l4),
    xlab = "Easting", ylab = "Northing", 
    scales = list(draw = TRUE),  at = breakpoints, 
    col.regions = terrain.colors(n = length(breakpoints)-1),
    par.settings=list(fontsize=list(text=20)))


###################################################
### code chunk number 15: CARBayesST.Rnw:669-673
###################################################
library("spdep")
W.nb <- poly2nb(GGHB.IG, row.names = SMR.av$IG)
W.list <- nb2listw(W.nb, style = "B")
W <- nb2mat(W.nb, style = "B")


###################################################
### code chunk number 16: CARBayesST.Rnw:683-689
###################################################
formula <- observed ~ offset(log(expected)) + jsa + price + pm10
model1 <- glm(formula = formula, family = "quasipoisson", 
    data = pollutionhealthdata)
resid.glm <- residuals(model1)
summary(model1)$coefficients
summary(model1)$dispersion


###################################################
### code chunk number 17: CARBayesST.Rnw:694-695
###################################################
moran.mc(x = resid.glm[1:271], listw = W.list, nsim = 10000)


###################################################
### code chunk number 18: CARBayesST.Rnw:820-824
###################################################
library("CARBayesdata")
library("sp")
data("GGHB.IG")
data("salesdata")


###################################################
### code chunk number 19: CARBayesST.Rnw:830-834
###################################################
salesdata$salesprop <- salesdata$sales / salesdata$stock
boxplot(salesdata$salesprop ~ salesdata$year, range = 0, xlab = "Year", 
    ylab = "Property sales rate", 
    col = "darkseagreen", border = "navy")


###################################################
### code chunk number 20: CARBayesST.Rnw:847-865
###################################################
library("dplyr")
salesprop.av <- summarise(group_by(salesdata, IG), 
    salesprop.mean = mean(salesprop))
GGHB.IG@data$sales <- salesprop.av$salesprop.mean
l1 = list("SpatialPolygonsRescale", layout.north.arrow(), 
    offset = c(220000,647000), scale = 4000)
l2 = list("SpatialPolygonsRescale", layout.scale.bar(), 
    offset = c(225000,647000), scale = 10000, 
    fill = c("transparent","black"))
l3 = list("sp.text", c(225000,649000), "0")
l4 = list("sp.text", c(230000,649000), "5000 m")
breakpoints <- c(0, quantile(salesprop.av$salesprop.mean, 
    seq(0.1, 0.9, 0.1)), 0.1)
spplot(GGHB.IG, "sales", sp.layout=list(l1, l2, l3, l4),
    xlab="Easting", ylab="Northing", 
    scales=list(draw = TRUE), at=breakpoints,
    col.regions=terrain.colors(n=length(breakpoints)-1),
    par.settings=list(fontsize=list(text=20)))


###################################################
### code chunk number 21: CARBayesST.Rnw:881-884
###################################################
library("spdep")
W.nb <- poly2nb(GGHB.IG, row.names = salesprop.av$salesprop.mean)
W <- nb2mat(W.nb, style = "B")


