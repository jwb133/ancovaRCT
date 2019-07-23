library(sandwich)
library(xtable)

set.seed(7912356)

#function to perform simulation study for a given scenario
performSim <- function(nSim=100000) {
  est <- array(0, dim=nSim)
  stdSE <- array(0, dim=nSim)
  sandwichSE <- array(0, dim=nSim)
  stdCI <- array(0, dim=nSim)
  sandwichCI <- array(0, dim=nSim)

  for (i in 1:nSim) {
    simDataReturn <- generateData()
    simData <- simDataReturn$data
    n <- dim(simData)[1]

    #model based analysis
    mod <- lm(y~a+w, data=simData)
    est[i] <- coef(mod)[2]
    stdSE[i] <- vcov(mod)[2,2]^0.5
    stdCI[i] <- 1*(((mod$coefficients[2]-qt(0.975,df=mod$df.residual)*summary(mod)$coefficients[2,2])<simDataReturn$trueVal) &
                     ((mod$coefficients[2]+qt(0.975,df=mod$df.residual)*summary(mod)$coefficients[2,2])>simDataReturn$trueVal))

    #robust sandwich variance and CI
    sandwichCov <- vcovHC(mod, type="HC3")
    sandwichSE[i] <- sandwichCov[2,2]^0.5
    sandwichCI[i] <- 1*(((mod$coefficients[2]-qt(0.975,df=mod$df.residual)*sandwichCov[2,2]^0.5)<simDataReturn$trueVal) &
                          ((mod$coefficients[2]+qt(0.975,df=mod$df.residual)*sandwichCov[2,2]^0.5)>simDataReturn$trueVal))

  } 
  #calculate ratio of mean variance estimates to empirical variance of estimated treatment effects
  #and coverage of CIs
  c(simDataReturn$varstarAncova, var(sqrt(n)*est), simDataReturn$nvarM, n*mean(stdSE^2), n*mean(sandwichSE^2),
      100*mean(stdCI), 100*mean(sandwichCI))
}

resTable <- array(0, dim=c(5,7))

#define five data generating mechanisms and for each call performSim(), saving results in resTable

#conditional mean correct, residual variance not dependent on A
generateData <- function() {
  n <- 500
  #higher probability of randomisation to experimental (a=1)
  pi <- 2/3
  w <- runif(n)
  a <- 1*(runif(n)<pi)
  y <- a+w+rnorm(n, sd=w^0.5)
  #var(Y-betaw*W|A)=var(eps)=E(var(eps|W))+var(E(eps|W))=E(W)=0.5
  list(data=data.frame(w=w,a=a,y=y), trueVal=1, varstarAncova=(0.5/pi+0.5/(1-pi)), nvarM=(0.5/(1-pi)+0.5/pi))
}
resTable[1,] <- performSim()

#conditional mean correct, conditional variance higher in experimental
generateData <- function() {
  n <- 500
  pi <- 2/3
  w <- runif(n)
  a <- 1*(runif(n)<pi)
  y <- a+w+rnorm(n, sd=(1+a)^0.5)
  #var(Y-betaw*W|A)=var(eps|A)=1+A
  list(data=data.frame(w=w,a=a,y=y), trueVal=1, varstarAncova=(2/pi+1/(1-pi)), nvarM=(2/(1-pi)+1/pi))
}
resTable[2,] <- performSim()

#conditional mean correct, conditional variance lower in experimental
generateData <- function() {
  n <- 500
  pi <- 2/3
  w <- runif(n)
  a <- 1*(runif(n)<pi)
  y <- a+w+rnorm(n, sd=(2-a)^0.5)
  #var(Y-betaw*W|A)=var(eps|A)=2-A
  list(data=data.frame(w=w,a=a,y=y), trueVal=1, varstarAncova=(1/pi+2/(1-pi)), nvarM=(1/(1-pi)+2/pi))
}
resTable[3,] <- performSim()

#conditional mean incorrect, but no interaction between a and w, and conditional variance not dependent on A
generateData <- function() {
  n <- 500
  pi <- 2/3
  w <- runif(n)
  a <- 1*(runif(n)<pi)
  y <- a+w+w^2+rnorm(n, sd=w^0.5)
  #some algebra with properties of the uniform distribution on (0,1) shows betaw_underbar=2 here
  #further algebra shows the Var(y-beta_w_underbar|A) is
  resVar <- 1/2 + (1/5 - 1/2 + 1/3 - (1/3-1/2)^2)

  list(data=data.frame(w=w,a=a,y=y), trueVal=1, varstarAncova=(resVar/pi+resVar/(1-pi)), nvarM=(resVar/(1-pi)+resVar/pi))
}
resTable[4,] <- performSim()

#conditional mean incorrect, interaction between x and z, residual variance same in both groups
generateData <- function() {
  n <- 500
  pi <- 2/3
  w <- runif(n)
  a <- 1*(runif(n)<pi)
  y <- a+4*a*w+rnorm(n, sd=1)
  #betaw_underbar is weighted average of 4 (wgt 2/3) and 0 (wgt 1/3) = 8/3
  resVar1 <- ((4/3)^2)*(1/12) + 1
  resVar0 <- ((8/3)^2)*(1/12) + 1

  list(data=data.frame(w=w,a=a,y=y), trueVal=3, varstarAncova=(resVar1/pi+resVar0/(1-pi)), nvarM=(resVar1/(1-pi)+resVar0/pi))
}
resTable[5,] <- performSim()

modelTable <- resTable[,c(1,2,3,4,6)]
xtable(modelTable)

sandwichTable <- resTable[,c(1,2,5,7)]
xtable(sandwichTable)
