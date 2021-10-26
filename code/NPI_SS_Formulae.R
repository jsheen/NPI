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
### Inputs: ###
## Rt1 is the effective reproductive number in the intervention group
## Rt0 is the effective reproductive number in the control group
## k1 is the overdispersion parameter for transmission in the intervention group
## k0 is the overdispersion parameter for transmission in the control group
## n is the population size for each cluster (assumed constant)
## EIt is the expected (mean across clusters) proportion of infectious individuals at time t, E[I_t].
## VarIt is the variance of the proportion of infectious individuals at time t, Var(I_t), across clusters. (This can be left as NULL)
## power, alpha, and N: specify 2 and leave 1 as NULL to solve for that one. 
#### Note that N (or the result if N is null) is the number of clusters in EACH arm

### Simplifications: ###
## Uses E[1/It] = 1/E[It], which underestimates variance
## Exactly one generation has passed between I_t and I_t+1
## Assumes full testing: everyone in the cluster is tested
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
## EIt is the expected (mean across clusters) proportion of infectious individuals at time t, E[I_t].
## power, alpha, and N: specify 2 and leave 1 as NULL to solve for that one. 
#### Note that N (or the result if N is null) is the number of clusters in EACH arm

### Simplifications: ###
## Disregards sampling variability in I_t
## Uses E[1/It] = 1/E[It], which underestimates variance
## Exactly one generation has passed between I_t and I_t+1
SScalc.samp <- function(Rt1, Rt0, k1, k0, m, n, EIt,
                   pow=.8, alpha=.05, N=NULL) {
  if (is.null(pow) + is.null(alpha) + is.null(N) != 1) {
    return("Error: Need exactly one of power, alpha, and N to be NULL.")
  } else {
    sd1 <- sqrt((Rt1/m)*(1/EIt-Rt1+(m-1)/n*(1+Rt1/k1)/EIt))
    sd0 <- sqrt((Rt0/m)*(1/EIt-Rt0+(m-1)/n*(1+Rt0/k0)/EIt))
    return(t.test.calcs(sd1=sd1, sd2=sd0, alpha=alpha, N=N, delta=Rt1-Rt0, pow=pow))
  }
}

