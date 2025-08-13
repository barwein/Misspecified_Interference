# Causal Inference with Misspecified Network Interference Structure
R package that implements the network-misspecification-robust (NMR) estimators. \
Theoretic background and full details are provided in the paper (https://arxiv.org/abs/2302.11322).

## Installation
For latest version please run

```{r}
install.packages("devtools")
devtools::install_github("barwein/Misspecified_Interference")
library(MisspecifiedInterference)
```

## Using the R package

The package contains one main function

-- `NMR_estimator` which estimate causal effects using the NMR HT and Hajek estimators.
