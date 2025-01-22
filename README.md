# NLBayes (R package)

A Bayesian Networks approach for inferring active Transcription Factors
using logic models of transcriptional regulation.

## Examples

The package includes example scripts:

- [nlbayes_example.R](examples/nlbayes_example.R) - Basic TF inference workflow using real data
- [tf_inference_simulation.R](examples/tf_inference_simulation.R) - Simulation-based example demonstrating:
  - Random network generation
  - Evidence simulation
  - OR-NOR model inference

## Installation Instructions

### Prerequisites

The package requires the GNU Scientific Library (GSL). Install it for your operating system:

#### Linux (Debian/Ubuntu)
```bash
# Install GSL
sudo apt install -y libgsl-dev

# Install devtools dependencies
sudo apt install -y libcurl4-openssl-dev libxml2-dev libssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

#### macOS
```bash
# Install GSL
brew install gsl
```

#### Windows
1. Download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
2. Install MSYS2 and open the MSYS2 terminal
3. Install GSL:
```bash
pacman -S mingw-w64-x86_64-gsl
```
4. Add the MSYS2 binary path (typically `C:\msys64\mingw64\bin`) to your system's PATH environment variable

### Package Installation

Install nlbayes directly from GitHub:
```r
# Install remotes if not already installed
install.packages("remotes")

# Install nlbayes
remotes::install_github("umbibio/nlbayes-r")
```

## Usage

Basic example of running TF inference:

```r
library(nlbayes)

# Create model
model <- ornor(
    network = network_list,  # List mapping TF IDs to target genes with regulation modes
    evidence = expression_data,  # Named vector of gene expression states
    n_graphs = 3,  # Number of parallel graphs for convergence assessment
    uniform_prior = FALSE,  # Whether to use uniform prior for theta parameter
    active_tfs = c("TF1", "TF2")  # Optional: TFs known to be active
)

# Fit the model
model <- fit(model,
    n_samples = 2000,  # Number of MCMC samples
    gelman_rubin = 1.1,  # Convergence criterion
    burnin = TRUE  # Whether to perform burn-in phase
)

# Get results
results <- get_results(model)
print(results)
```

The input data should be formatted as follows:

- `network_list`: A named list where names are TF IDs and values are named numeric vectors mapping target gene names to regulation modes (-1 for repression, 1 for activation)
- `expression_data`: A named numeric vector of gene expression states (-1 for down-regulated, 1 for up-regulated)

For more detailed examples, check the example scripts in the `examples` directory.
