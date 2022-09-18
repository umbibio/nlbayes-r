# Noisy Logic Bayesian models for active TF inference

## Install on Linux

1. Install system dependencies

Devtools dependencies
```
$ sudo apt install -y libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```
GNU gsl library (nlbayes dependency)
```
$ sudo apt install -y libgsl-dev
```
2. Install R dependencies
```
$ R -q
> install.packages(c('Rcpp', 'RcppProgress', 'devtools'))
```
3. Install the NLBayes package
```
> devtools::install_github('umbibio/nlbayes-rcran')
```
