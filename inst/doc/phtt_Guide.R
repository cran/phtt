### R code from vignette source 'phtt_Guide.Snw'

###################################################
### code chunk number 1: phtt_Guide.Snw:147-176
###################################################
## Install package
## install.packages("phtt", repos="http://R-Forge.R-project.org")
## Load package
library("phtt")
## vignette()
## vignette("phtt_Guide")
## Load Data
data("Cigar")
N <- 46
T <- 30
## Dependent variable:
## Cigarette-Sales per Capita
l.Consumption    <- log(matrix(Cigar$sales, T, N))
## Independent variables:
## Consumer Price Index
cpi              <- matrix(Cigar$cpi, T, N)
## Real Price per Pack of Cigarettes 
l.Price          <- log(matrix(Cigar$price, T, N)/cpi)
## Real Disposable Income per Capita  
l.Income         <- log(matrix(Cigar$ndi, T, N)/cpi)

pdf("Cigar.pdf")
scl <- 1.6
par(mfrow=c(1,3), mar=c(6, 5, 5, 2.1))
matplot(l.Consumption,type="l", main="Log's of\nCigar-Consumption",ylab="",xlab="Time",cex.main=scl,cex.lab=scl,cex.axis=scl)
matplot(l.Price, type="l",      main="Log's of\nreal Prices",    ylab="",xlab="Time",cex.main=scl,cex.lab=scl, cex.axis=scl)
matplot(l.Income, type="l",     main="Log's of\nreal Income",    ylab="",xlab="Time",cex.main=scl,cex.lab=scl, cex.axis=scl)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
graphics.off()


###################################################
### code chunk number 2: phtt_Guide.Snw:291-292
###################################################
args(KSS)


###################################################
### code chunk number 3: phtt_Guide.Snw:316-324
###################################################
library("phtt")
data("Cigar")
N <- 46
T <- 30
l.Consumption   <- log(matrix(Cigar$sales, T, N))
cpi             <- matrix(Cigar$cpi,       T, N)
l.Price         <- log(matrix(Cigar$price, T, N)/cpi)
l.Income        <- log(matrix(Cigar$ndi,   T, N)/cpi)


###################################################
### code chunk number 4: phtt_Guide.Snw:329-331
###################################################
Cigar.KSS <- KSS(formula = l.Consumption ~ l.Price + l.Income) 
(Cigar.KSS.summary <- summary(Cigar.KSS))


###################################################
### code chunk number 5: phtt_Guide.Snw:335-340
###################################################
## Figure 2:
pdf("KSSM1.pdf")
scl <- 1
plot(Cigar.KSS.summary,cex.main=scl,cex.lab=scl,cex.axis=scl)
graphics.off()


###################################################
### code chunk number 6: phtt_Guide.Snw:344-345 (eval = FALSE)
###################################################
## plot(Cigar.KSS.summary)


###################################################
### code chunk number 7: phtt_Guide.Snw:437-438
###################################################
args(OptDim)


###################################################
### code chunk number 8: phtt_Guide.Snw:448-449
###################################################
OptDim(Obj = l.Consumption, criteria = "PC1")


###################################################
### code chunk number 9: phtt_Guide.Snw:454-456
###################################################
(OptDim.obj <- OptDim(Obj = l.Consumption, criteria = c("PC3",  "ER",  
                      "GR", "IPC1", "IPC2", "IPC3"), standardize = TRUE))


###################################################
### code chunk number 10: phtt_Guide.Snw:461-464
###################################################
pdf("OptDimv.pdf")
plot(OptDim.obj)
graphics.off()


###################################################
### code chunk number 11: phtt_Guide.Snw:468-469 (eval = FALSE)
###################################################
## plot(OptDim.obj) 


###################################################
### code chunk number 12: phtt_Guide.Snw:485-486 (eval = FALSE)
###################################################
## KSS(formula = l.Consumption ~ -1 + l.Price + l.Income, consult.dim = TRUE) 


###################################################
### code chunk number 13: phtt_Guide.Snw:676-677
###################################################
args(Eup)


###################################################
### code chunk number 14: phtt_Guide.Snw:692-695
###################################################
d.l.Consumption  <- diff(l.Consumption)
d.l.Price        <- diff(l.Price)
d.l.Income       <- diff(l.Income)


###################################################
### code chunk number 15: phtt_Guide.Snw:699-701
###################################################
(Cigar.Eup <- Eup(d.l.Consumption ~  -1 + d.l.Price + d.l.Income, 
                  dim.criterion = "PC3"))


###################################################
### code chunk number 16: phtt_Guide.Snw:717-718
###################################################
summary(Cigar.Eup)


###################################################
### code chunk number 17: phtt_Guide.Snw:724-727
###################################################
pdf("EupPlot.pdf")
plot(summary(Cigar.Eup))
graphics.off()


###################################################
### code chunk number 18: phtt_Guide.Snw:731-732 (eval = FALSE)
###################################################
## plot(summary(Cigar.Eup))


###################################################
### code chunk number 19: phtt_Guide.Snw:797-800
###################################################
Cigar2.KSS <- KSS(formula = l.Consumption ~ l.Price + l.Income,
                  additive.effects = "individual") 
Cigar2.KSS.summary <- summary(Cigar2.KSS)


###################################################
### code chunk number 20: phtt_Guide.Snw:803-806 (eval = FALSE)
###################################################
## Cigar2.KSS <- KSS(formula = l.Consumption ~ l.Price + l.Income,
##                   additive.effects = "individual") 
## (Cigar2.KSS.summary <- summary(Cigar2.KSS))


###################################################
### code chunk number 21: phtt_Guide.Snw:838-839 (eval = FALSE)
###################################################
## plot(Cigar2.KSS.summary)


###################################################
### code chunk number 22: phtt_Guide.Snw:842-848
###################################################
pdf("KSSM2.pdf")
scl <- 1.75
par(mar=c(6, 5, 5, 2.1))
plot(Cigar2.KSS.summary,cex.main=1.25,cex.lab=scl,cex.axis=scl)
par(mar=c(5.1, 4.1, 4.1, 2.1))
graphics.off()


###################################################
### code chunk number 23: phtt_Guide.Snw:892-893 (eval = FALSE)
###################################################
## checkSpecif(obj1, obj2, level = 0.05)


###################################################
### code chunk number 24: phtt_Guide.Snw:907-912 (eval = FALSE)
###################################################
## twoways.obj     <- Eup(d.l.Consumption ~  -1 + d.l.Price + d.l.Income, 
##                        factor.dim = 0, additive.effects = "twoways")
## not.twoways.obj <- Eup(d.l.Consumption ~  -1 + d.l.Price + d.l.Income, 
##                        factor.dim = 2, additive.effects = "none")  
## checkSpecif(obj1 = twoways.obj, obj2 = not.twoways.obj, level = 0.01)


###################################################
### code chunk number 25: phtt_Guide.Snw:953-956
###################################################
Eup.obj <- Eup(d.l.Consumption ~  -1 + d.l.Price + d.l.Income, 
               additive.effects = "twoways")
checkSpecif(Eup.obj, level = 0.01)


###################################################
### code chunk number 26: phtt_Guide.Snw:960-963
###################################################
KSS.obj <- KSS(l.Consumption ~  -1 + l.Price + l.Income, 
               additive.effects = "twoways")
checkSpecif(KSS.obj, level = 0.01)


###################################################
### code chunk number 27: phtt_Guide.Snw:1010-1021 (eval = FALSE)
###################################################
## coef(Cigar2.KSS)$Var.shares.of.loadings.param[1]
## coef(Cigar2.KSS)$Var.shares.of.loadings.param[2]
## coef(Cigar2.KSS)$Var.shares.of.loadings.param[3]
## coef(Cigar2.KSS)$Var.shares.of.loadings.param[4]
## coef(Cigar2.KSS)$Var.shares.of.loadings.param[5]
## 
## coef(Cigar2.KSS)$Common.factors[,2]
## lambda_i1 <- coef(Cigar2.KSS)$Ind.loadings.param[,1]
## order(abs(lambda_i1),decreasing=TRUE)[1]
## coef(Cigar2.KSS)$Ind.loadings.param[7,1]
## round(range(coef(Cigar2.KSS)$Ind.loadings.param[-7,1]), digits=2)


###################################################
### code chunk number 28: phtt_Guide.Snw:1024-1031
###################################################
pdf("Factor2.pdf")
scl <- 1
par(mfrow=c(1,2))
matplot(coef(Cigar2.KSS)$Common.factors[,1]%*%t(coef(Cigar2.KSS)$Ind.loadings.param[,1]),type="l",main="Variance of time-varying indiv. effects\n in direction of the 1. common factor",xlab="Time",ylab="",cex.main=scl,cex.lab=scl,cex.axis=scl)
matplot(coef(Cigar2.KSS)$Common.factors[,2]%*%t(coef(Cigar2.KSS)$Ind.loadings.param[,2]),type="l",main="Variance of time-varying indiv. effects\n in direction of the 2. common factor",xlab="Time",ylab="",cex.main=scl,cex.lab=scl,cex.axis=scl)
par(mfrow=c(1,1))
graphics.off()


