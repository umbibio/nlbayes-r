# NLBayes (R package)

Noisy Logic Bayesian models for active TF inference

To learn how to use this tool, check out the [example script](examples/nlbayes_example.R).


## Install on Linux

1. Install system dependencies

Devtools dependencies
```bash
sudo apt install -y libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```
GNU gsl library (nlbayes dependency)
```bash
sudo apt install -y libgsl-dev
```
2. Install R dependencies
```bash
R -q -e "install.packages(c('rjson', 'Rcpp', 'RcppProgress', 'devtools'))"
```
3. Install the NLBayes package
```bash
R -q -e "devtools::install_github('umbibio/nlbayes-rcran')"
```
