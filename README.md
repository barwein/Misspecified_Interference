# Causal inference with misspecified interference structure
R package that implements the network-misspecification-robust (NMR) estimator and the sensitivity analysis for network misspecification. \
Theoretic background and full details are provided in the paper ``Causal inference with misspecified interference structure".

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
2. `SensitivityAnalysis` which run the sensitivity analysis for network misspecification (currently only supports contamination in CRT).
