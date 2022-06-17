# PowerAnalysis R Package

Helloword

#### Power/Precision Analysis to Determine Minimum Sample Size and Effect

The PowerAnalysis R package contains functions to perform power analysis and/or
precision analysis for sample means of continuous and binary metrics and sample
quantiles of continuous metrics, assuming asymptotic normality. One can
determine the minimum required sample size (MRSS), the minimum required duration
(MRD), or the minimum detectable effect (MDE) and account for design effects
due to arm imbalance, analysis of covariance (ANCOVA), and cluster
randomization. There are also helper functions to forecast sample size growth
over time to convert between sample size growth and duration.

#### Installation

```r
install.packages("devtools")
library(devtools)
devtools::install_github("google/PowerAnalysis")
library(PowerAnalysis)
```

#### Disclaimer

This is not an officially supported Google product.
