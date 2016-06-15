### R code from vignette source 'CARBayesST.Rnw'

###################################################
### code chunk number 1: CARBayesST.Rnw:95-96
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 2: CARBayesST.Rnw:357-363
###################################################
x.easting <- 1:20
x.northing <- 1:20
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)
N <- 20
N.all <- N * K


###################################################
### code chunk number 3: CARBayesST.Rnw:368-377
###################################################
W <-array(0, c(K,K))
    for(i in 1:K)
    {
        for(j in 1:K)
        {
        temp <- (Grid[i,1] - Grid[j,1])^2 + (Grid[i,2] - Grid[j,2])^2
            if(temp==1)  W[i,j] <- 1 
        }    
    }


###################################################
### code chunk number 4: CARBayesST.Rnw:382-390
###################################################
D <-array(0, c(N,N))
    for(i in 1:N)
    {
        for(j in 1:N)
        {
            if(abs((i-j))==1)  D[i,j] <- 1 
        }    
    }


###################################################
### code chunk number 5: CARBayesST.Rnw:395-396
###################################################
Q.W <- 0.99 * (diag(apply(W, 2, sum)) - W) + 0.01 * diag(rep(1,K))


###################################################
### code chunk number 6: CARBayesST.Rnw:401-405
###################################################
Q.W.inv <- solve(Q.W)
library(MASS)
phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.01 * Q.W.inv))
phi.long <- rep(phi, N)


###################################################
### code chunk number 7: CARBayesST.Rnw:411-415
###################################################
Q.D <- 0.99 * (diag(apply(D, 2, sum)) - D) + 0.01 * diag(rep(1,N))
Q.D.inv <- solve(Q.D)
delta <- mvrnorm(n=1, mu=rep(0,N), Sigma=(0.01 * Q.D.inv))
delta.long <- kronecker(delta, rep(1,K))


###################################################
### code chunk number 8: CARBayesST.Rnw:420-422
###################################################
x <- rnorm(n=N.all, mean=0, sd=1)
gamma <- rnorm(n=N.all, mean=0, sd=0.1)


###################################################
### code chunk number 9: CARBayesST.Rnw:427-433
###################################################
beta1 <- -0.1
beta2 <- 0.1
n <- rep(200, N.all)
LP <- beta1 + beta2 * x + phi.long +  delta.long + gamma
theta.true <- exp(LP) / (1 + exp(LP))
Y <- rbinom(n=N.all, size=n, prob=theta.true)


###################################################
### code chunk number 10: CARBayesST.Rnw:521-525
###################################################
library(CARBayesdata)
library(sp)
data(GGHB.IG)
data(pollutionhealthdata)


###################################################
### code chunk number 11: CARBayesST.Rnw:531-532
###################################################
head(pollutionhealthdata)


###################################################
### code chunk number 12: CARBayesST.Rnw:538-542
###################################################
pollutionhealthdata$SMR <- pollutionhealthdata$observed / pollutionhealthdata$expected
pollutionhealthdata$logSMR <- log(pollutionhealthdata$SMR)
pairs(pollutionhealthdata[ ,c(9, 5:7)], pch=19, cex=0.5, lower.panel=NULL,
      labels=c("ln(SMR)", "PM10", "JSA", "Price (*100,000)"))


###################################################
### code chunk number 13: CARBayesST.Rnw:555-558
###################################################
library(dplyr)
SMR.av <- summarise(group_by(pollutionhealthdata,IG), SMR.mean=mean(SMR))
GGHB.IG@data$SMR <- SMR.av$SMR.mean


###################################################
### code chunk number 14: CARBayesST.Rnw:564-575
###################################################
l1 = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000,
    647000), scale = 4000)
l2 = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000,
    647000), scale = 10000, fill=c("transparent","black"))
l3 = list("sp.text", c(225000,649000), "0")
l4 = list("sp.text", c(230000,649000), "5000 m")
breakpoints <- seq(min(SMR.av$SMR.mean)-0.1, max(SMR.av$SMR.mean)+0.1, length.out=11)
spplot(GGHB.IG, "SMR", sp.layout=list(l1, l2, l3, l4),
       col="transparent", xlab="Easting", ylab="Northing", 
       scales=list(draw = TRUE), at=breakpoints,
       col.regions=terrain.colors(n=length(breakpoints)-1))


###################################################
### code chunk number 15: CARBayesST.Rnw:588-592
###################################################
library(spdep)
W.nb <- poly2nb(GGHB.IG, row.names = SMR.av$IG)
W.list <- nb2listw(W.nb, style="B")
W <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 16: CARBayesST.Rnw:602-607
###################################################
formula <- observed~offset(log(expected)) + jsa + price + pm10
model1 <- glm(formula=formula, family="quasipoisson", data=pollutionhealthdata)
resid.glm <- residuals(model1)
summary(model1)$coefficients
summary(model1)$dispersion


###################################################
### code chunk number 17: CARBayesST.Rnw:612-613
###################################################
moran.mc(x=resid.glm[1:271], listw = W.list, nsim = 10000)


###################################################
### code chunk number 18: CARBayesST.Rnw:738-742
###################################################
library(CARBayesdata)
library(sp)
data(GGHB.IG)
data(salesdata)


###################################################
### code chunk number 19: CARBayesST.Rnw:748-751
###################################################
salesdata$salesprop <- salesdata$sales / salesdata$stock
boxplot(salesdata$salesprop~salesdata$year, range=0, xlab="Year", ylab="Property sales 
        as a proportion of total properties", col="darkseagreen", border="navy")


###################################################
### code chunk number 20: CARBayesST.Rnw:764-778
###################################################
library(dplyr)
salesprop.av <- summarise(group_by(salesdata,IG), salesprop.mean=mean(salesprop))
GGHB.IG@data$sales <- salesprop.av$salesprop.mean
l1 = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(220000,647000), 
          scale = 4000)
l2 = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(225000,647000), 
          scale = 10000, fill=c("transparent","black"))
l3 = list("sp.text", c(225000,649000), "0")
l4 = list("sp.text", c(230000,649000), "5000 m")
breakpoints <- c(0, quantile(salesprop.av$salesprop.mean, seq(0.1,0.9,0.1)), 0.1)
spplot(GGHB.IG, "sales", sp.layout=list(l1, l2, l3, l4),
       col="transparent", xlab="Easting", ylab="Northing", 
       scales=list(draw = TRUE), at=breakpoints,
       col.regions=terrain.colors(n=length(breakpoints)-1))


###################################################
### code chunk number 21: CARBayesST.Rnw:794-797
###################################################
library(spdep)
W.nb <- poly2nb(GGHB.IG, row.names = salesprop.av$salesprop.mean)
W <- nb2mat(W.nb, style="B")


