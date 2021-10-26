### Requires tidyverse package and NPI_SS_Formulae.R commands
require(tidyverse)
source(file="NPI_SS_Formulae.R")

## Proposed Approach Section Approximations:
### Full testing:
SScalc.full(Rt1=1.2*.6, Rt0=1.2, k1=0.4, k0=0.4, n=570, EIt=0.005)

### Testing 100 individuals per cluster:
SScalc.samp(Rt1=1.2*.6, Rt0=1.2, k1=0.4, k0=0.4, m=100, n=570, EIt=0.005)

### Testing 50 individuals per cluster:
SScalc.samp(Rt1=1.2*.6, Rt0=1.2, k1=0.4, k0=0.4, m=50, n=570, EIt=0.005)

### Sensitivity to k parameter:
SScalc.full(Rt1=1.2*.6, Rt0=1.2, k1=0.1, k0=0.1, n=570, EIt=0.005)
SScalc.samp(Rt1=1.2*.6, Rt0=1.2, k1=0.1, k0=0.1, m=100, n=570, EIt=0.005)
SScalc.samp(Rt1=1.2*.6, Rt0=1.2, k1=0.1, k0=0.1, m=50, n=570, EIt=0.005)

### Sensitivity to R parameter:
SScalc.full(Rt1=1.5*.6, Rt0=1.5, k1=0.4, k0=0.4, n=570, EIt=0.005)
SScalc.samp(Rt1=1.5*.6, Rt0=1.5, k1=0.4, k0=0.4, m=100, n=570, EIt=0.005)
SScalc.samp(Rt1=1.5*.6, Rt0=1.5, k1=0.4, k0=0.4, m=50, n=570, EIt=0.005)
