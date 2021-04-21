require(tidyverse)

### Function to calculate power, effect size, or sample size per arm in two-sided two-sample t-test 
#### with equal arm sizes and unequal variances:
t.test.calcs <- function(sd1, sd2, alpha=0.05, N=NULL, delta=NULL, pow=NULL) {
  if (is.null(N)+is.null(delta)+is.null(pow)+is.null(alpha) != 1) {
    print("Error: Need exactly one of N, delta, power, or alpha to be null")
    return(NA)
  } else if (is.null(N)) {
    approx <- ((sd1^2+sd2^2)*(qnorm(1-alpha/2)+qnorm(pow))^2)/(delta^2)
    if (approx < 1) {
      return(1)
    } else {
      return(uniroot(f=function(x) ((sd1^2+sd2^2)*(qt(p=1-alpha/2, df=2*x-2)+qt(p=pow, df=2*x-2))^2)/(delta^2)-x,
                   lower=approx, upper=approx*5)$root)
    }
  } else if(is.null(pow)) {
    return(pow <- pt(q=sqrt((N*delta^2)/(sd1^2+sd2^2))-qt(p=1-alpha/2, df=2*N-2), df=2*N-2))
  } else if(is.null(delta)) {
    return(sqrt(((sd1^2+sd2^2)*(qt(p=1-alpha/2, df=2*N-2)+qt(p=pow, df=2*N-2))^2)/N))
  } else {
    return(2*(1-pt(q=sqrt((N*delta^2)/(sd1^2+sd2^2))-qt(p=pow, df=2*N-2), df=2*N-2)))
  }
}

### Full Testing: ###
SScalc.full <- function(Rt1, Rt0, k1, k0, n, EIt, VarIt=NULL,
                        pow=.8, alpha=.05, N=NULL) {
  if (is.null(pow) + is.null(alpha) + is.null(N) != 1) {
    return("Error: Need exactly one of power, alpha, and N to be NULL.")
  } else {
    sd1 <- sqrt(Rt1*(1+Rt1/k1)/(n*EIt))
    sd0 <- sqrt(Rt0*(1+Rt0/k0)/(n*EIt))
    calc1 <- t.test.calcs(sd1=sd1, sd2=sd0, alpha=alpha, N=N, delta=Rt1-Rt0, pow=pow)
    if (is.null(VarIt)) {
      return(calc1)
    } else {
      sd1 <- sqrt(Rt1*(1+Rt1/k1)*(1/(n*EIt)+VarIt/(n*EIt^3)))
      sd0 <- sqrt(Rt0*(1+Rt0/k0)*(1/(n*EIt)+VarIt/(n*EIt^3)))
      calc2 <- t.test.calcs(sd1=sd1, sd2=sd0, alpha=alpha, N=N, delta=Rt1-Rt0, pow=pow)
      return(c(Lower=calc1,Upper=calc2))
    }
  }
}

### Inputs: ###
## Rt1 is the effective reproductive number in the intervention group
## Rt0 is the effective reproductive number in the control group
## k1 is the overdispersion parameter for transmission in the intervention group
## k0 is the overdispersion parameter for transmission in the control group
## m is the number of individuals sampled in each cluster at time t+1
## n is the population size for each cluster (assumed constant)
## EIt is the expected number of infectious individuals at time t, E[I_t].
## power, alpha, and N: specify 2 and leave 1 as NULL to solve for that one. 
#### Note that N (or the result if N is null) is the number of clusters in EACH arm

### Simplifications: ###
## Disregards sampling variability in I_t
## Uses E[1/It] = E[It], which reduces variance (lower bound on required sample size)
## Exactly one generation has passed between I_t and I_t+1
SScalc.LB <- function(Rt1, Rt0, k1, k0, m, n, EIt,
                   pow=.8, alpha=.05, N=NULL) {
  if (is.null(pow) + is.null(alpha) + is.null(N) != 1) {
    return("Error: Need exactly one of power, alpha, and N to be NULL.")
  } else {
    sd1 <- sqrt((Rt1/m)*(1/EIt-Rt1+(m-1)/n*(1+Rt1/k1)/EIt))
    sd0 <- sqrt((Rt0/m)*(1/EIt-Rt0+(m-1)/n*(1+Rt0/k0)/EIt))
    return(t.test.calcs(sd1=sd1, sd2=sd0, alpha=alpha, N=N, delta=Rt1-Rt0, pow=pow))
  }
}

### Second version, with a log-ratio test statistic,
### incorporating sampling at both time points but with other approximations:
SScalc.log <- function(Rt1, Rt0, k1, k0, m, m0, n, EIt,
                      pow=.8, alpha=.05, N=NULL) {
  if (is.null(pow) + is.null(alpha) + is.null(N) != 1) {
    return("Error: Need exactly one of power, alpha, and N1 to be NULL.")
  } else {
    sd.simple.1 <- sqrt(1/(m*Rt1*EIt)-1/m+1/(m0*EIt)-1/m0)
    sd.simple.0 <- sqrt(1/(m*Rt0*EIt)-1/m+1/(m0*EIt)-1/m0)
    res1 <- t.test.calcs(sd1=sd.simple.1, sd2=sd.simple.0, alpha=alpha, N=N, delta=log(Rt1/Rt0), pow=pow)
    sd1a <- (1/m-1/m^2)*(1/(Rt1*EIt)+(1+Rt1/k1)/(n*(Rt1*EIt)^2))
    sd1b <- (1/m^2)*(1/(Rt1*EIt)^2+3*(1+Rt1/k1)/(n*(Rt1*EIt)^3))-1/m
    sd1c <- 1/(m0*EIt)+(1/(m0^2*EIt))*(1/EIt-1)-1/m0
    sd1d <- (1+Rt1/k1)/(n*Rt1*EIt)
    sd1 <- sqrt(sd1a+sd1b+sd1c+sd1d)
    sd0a <- (1/m-1/m^2)*(1/(Rt0*EIt)+(1+Rt0/k0)/(n*(Rt0*EIt)^2))
    sd0b <- (1/m^2)*(1/(Rt0*EIt)^2+3*(1+Rt0/k0)/(n*(Rt0*EIt)^3))-1/m
    sd0c <- 1/(m0*EIt)+(1/(m0^2*EIt))*(1/EIt-1)-1/m0
    sd0d <- (1+Rt0/k0)/(n*Rt0*EIt)
    sd0 <- sqrt(sd0a+sd0b+sd0c+sd0d)
    res2 <- t.test.calcs(sd1=sd1, sd2=sd0, alpha=alpha, N=N, delta=log(Rt1/Rt0), pow=pow)
    return(list(Full=res2, Simple=res1))
  }
}

SScalc.log.opt <- function(Rt1, Rt0, k1, k0, m.total, n, EIt,
                          pow=.8, alpha=.05, type="Full") {
 DF <- data.frame(m=seq(1,m.total-1,by=1), m0=seq(m.total-1,1,by=-1))
 DF$N.full <- apply(DF, 1, function(x) SScalc.log(Rt1, Rt0, k1, k0, x[1], x[2],
                                            n, EIt, pow, alpha, N=NULL)$Full)
 DF$N.simple <- apply(DF, 1, function(x) SScalc.log(Rt1, Rt0, k1, k0, x[1], x[2],
                                                    n, EIt, pow, alpha, N=NULL)$Simple)
 return(rbind(Full=DF[DF$N.full==min(DF$N.full),],Simple=DF[DF$N.simple==min(DF$N.simple),]))
}

## E.g., for R_t under control of 2.0 and 40% reduction in R_t under intervention, 
### where k=0.5 in both intervention and control,
### where each cluster has a population of 10,000, we sample 100 from each cluster,
### and the infection rate at time t is, on average, 1%,
### we can find a lower bound for the number of clusters in each arm needed for 80% power and alpha=0.05:
SScalc.LB(Rt1=1.2, Rt0=2.0, k1=0.5, k0=0.5,
       m=100, n=10000, EIt=0.01,
       pow=.8, alpha=.05, N=NULL)
## ~42  clusters per arm.
SScalc.log(Rt1=1.2, Rt0=2.0, k1=0.5, k0=0.5,
          m=100, m0=100, n=10000, EIt=0.01,
          pow=.8, alpha=.05, N=NULL) ## Gives 100 clusters/arm
### For fixed m+m0, when Rt > 0, the best sampling is to sample more at t+1 than at t:
SScalc.log.opt(Rt1=1.2, Rt0=2.0, k1=0.5, k0=0.5,
              m.total=200, n=10000, EIt=0.01,
              pow=.8, alpha=.05)
## Use m=90, m0=110, with N=99 clusters/arm
##If we could test everyone and the sd of the infection rate at time t was .005:
SScalc.full(Rt1=1.2, Rt0=2.0, k1=0.5, k0=0.5,
            n=10000, EIt=0.01, VarIt=0.005^2,
            pow=.8, alpha=.05, N=NULL) ## 3-4 clusters/arm (too small for asymptotic t-test though!)

## Lower Rt:
SScalc.LB(Rt1=0.6, Rt0=1.0, k1=0.5, k0=0.5,
          m=100, n=10000, EIt=0.01,
          pow=.8, alpha=.05, N=NULL) ## 81 clusters/arm
SScalc.log(Rt1=0.6, Rt0=1.0, k1=0.5, k0=0.5,
          m=100, m0=100, n=10000, EIt=0.01,
          pow=.8, alpha=.05, N=NULL) ## 141 clusters/arm
SScalc.log.opt(Rt1=0.6, Rt0=1.0, k1=0.5, k0=0.5,
          m.total=200, n=10000, EIt=0.01,
          pow=.8, alpha=.05) ## Optimum: m=107, m0=93, N=140 clusters/arm
## With Rt < 1, it switches and there's more gain from sampling more at time t+1

## In-Between Rt:
SScalc.LB(Rt1=0.75, Rt0=1.25, k1=0.5, k0=0.5,
       m=100, n=10000, EIt=0.01,
       pow=.8, alpha=.05, N=NULL) ## 65 clusters/arm
SScalc.log(Rt1=0.75, Rt0=1.25, k1=0.5, k0=0.5,
          m=100, m0=100, n=10000, EIt=0.01,
          pow=.8, alpha=.05, N=NULL) ## 125 clusters/arm
SScalc.log.opt(Rt1=0.75, Rt0=1.25, k1=0.5, k0=0.5,
          m.total=200, n=10000, EIt=0.01,
          pow=.8, alpha=.05) ## m=102, m0=98, 124 clusters/arm
## With Rt straddling 1, ~best to evenly sample at t and t+1



### Comparison to Justin's Simulations:
Sim_Comps <- tibble(Rt0=rep(c(1.5,2),each=8), effect=rep(rep(c(.2,.4),each=4),2),
                    k=rep(c(.7,.7,.4,.4),4), popsize=rep(c(1000,10000),8))
Sim_Comps$full <- apply(Sim_Comps, MARGIN=1,
                        FUN=function(vals) SScalc.full(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"], 
                                                    k1=vals["k"], k0=vals["k"], n=vals["popsize"], 
                                                    EIt=.005, VarIt=NULL,
                                                    pow=.8, alpha=.05, N=NULL))
Sim_Comps$LB <- apply(Sim_Comps, MARGIN=1,
                      FUN=function(vals) SScalc.LB(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"], 
                                                     k1=vals["k"], k0=vals["k"], n=vals["popsize"], m=100,
                                                     EIt=.005,
                                                     pow=.8, alpha=.05, N=NULL))
Sim_Comps$log <- apply(Sim_Comps, MARGIN=1,
                      FUN=function(vals) SScalc.log(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"], 
                                                     k1=vals["k"], k0=vals["k"], n=vals["popsize"], m=100, m0=100,
                                                     EIt=.005,
                                                     pow=.8, alpha=.05, N=NULL)$Simple)

Sim_Comps$log2 <- apply(Sim_Comps, MARGIN=1,
                       FUN=function(vals) SScalc.log(Rt1=vals["Rt0"]*(1-vals["effect"]), Rt0=vals["Rt0"], 
                                                     k1=vals["k"], k0=vals["k"], n=vals["popsize"], m=100, m0=100,
                                                     EIt=.005,
                                                     pow=.8, alpha=.05, N=NULL)$Full)

#### Old Second version:
# SScalc.log <- function(Rt1, Rt0, k1, k0, m, m0, n, EIt,
#                        power=.8, alpha=.05, N=NULL) {
#   if (is.null(power) + is.null(alpha) + is.null(N) != 1) {
#     return("Error: Need exactly one of power, alpha, and N1 to be NULL.")
#   } else {
#     sd1 <- sqrt(1/(m*Rt1*EIt)-1/m+1/(m0*EIt)-1/m0)
#     sd0 <- sqrt(1/(m*Rt0*EIt)-1/m+1/(m0*EIt)-1/m0)
#     return(t.test.calcs(sd1=sd1, sd2=sd0, alpha=alpha, N=N, delta=log(Rt1/Rt0), power=power))
#   }
# }

SScalc.LB(Rt1=0.881, Rt0=1.5, k1=0.1, k0=0.1,
          m=1000, n=10000, EIt=0.005,
          pow=.8, alpha=.05, N=NULL)
