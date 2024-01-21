# Causal inference with misspecified interference structure
R package that implements the network-misspecification-robust (NMR) estimator and probabilistic bias analysis (PBA) of network misspecification. \
Theoretic background and full details are provided in the paper ``Causal inference with misspecified interference structure" (https://arxiv.org/abs/2302.11322).

## Installation
For latest version please run

```{r}
install.packages("devtools")
devtools::install_github("barwein/Misspecified_Interference")
library(MisspecifiedInterference)
```

## Using the R package

The package contains two functions

1. `NMR_estimator` which estimate causal effects using the NMR HT and Hajek estimators.
2. `PBA` which run the PBA for network misspecification with use choice of network deviation distribution and prior.
