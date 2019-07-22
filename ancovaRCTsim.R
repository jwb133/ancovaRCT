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
  c(mean(stdSE^2)/var(est), mean(sandwichSE^2)/var(est), 100*mean(stdCI), 100*mean(sandwichCI))
}

resTable <- array(0, dim=c(5,4))

#define five data generating mechanisms and for each call performSim(), saving results in resTable

#conditional mean correct, residual variance not dependent on A
generateData <- function() {
  n <- 500
  #higher probability of randomisation to experimental (a=1)
  pi <- 2/3
  w <- runif(n)
  a <- 1*(runif(n)<pi)
  y <- a+w+rnorm(n, sd=w^0.5)
  list(data=data.frame(w=w,a=a,y=y), trueVal=1)
}
resTable[1,] <- performSim()

#conditional mean correct, conditional variance higher in experimental
generateData <- function() {
  n <- 500
  pi <- 2/3
  w <- runif(n)
  a <- 1*(runif(n)<pi)
  y <- a+w+rnorm(n, sd=(1+a)^0.5)
  list(data=data.frame(w=w,a=a,y=y), trueVal=1)
}
resTable[2,] <- performSim()

#conditional mean correct, conditional variance lower in experimental
generateData <- function() {
  n <- 500
  pi <- 2/3
  w <- runif(n)
  a <- 1*(runif(n)<pi)
  y <- a+w+rnorm(n, sd=(2-a)^0.5)

  list(data=data.frame(w=w,a=a,y=y), trueVal=1)
}
resTable[3,] <- performSim()

#conditional mean incorrect, but no interaction between a and w, and conditional variance not dependent on A
generateData <- function() {
  n <- 500
  pi <- 2/3
  w <- runif(n)
  a <- 1*(runif(n)<pi)
  y <- a+w+w^2+rnorm(n, sd=w^0.5)

  list(data=data.frame(w=w,a=a,y=y), trueVal=1)
}
resTable[4,] <- performSim()

#conditional mean incorrect, interaction between x and z, residual variance same in both groups
generateData <- function() {
  n <- 500
  pi <- 2/3
  w <- runif(n)
  a <- 1*(runif(n)<pi)
  y <- a+4*a*w+rnorm(n, sd=1)

  list(data=data.frame(w=w,a=a,y=y), trueVal=3)
}
resTable[5,] <- performSim()

resTable

xtable(resTable)
