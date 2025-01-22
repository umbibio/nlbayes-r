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
# Install devtools if not already installed
install.packages("devtools")

# Install nlbayes (this will automatically install all required dependencies)
devtools::install_github('umbibio/nlbayes-r')
```

## Usage

Basic example of running TF inference:

```r
library(nlbayes)

# Load or generate your network and evidence
network <- list()  # List mapping TFs to their target genes
evidence <- c()    # Named numeric vector of differential expression evidence

# Create and run the model
inference.model <- ORNOR.inference(network, evidence, n.graphs = 5)
inference.model <- sample.posterior(inference.model, N = 2000, gr.level = 1.1, burnin = TRUE)

# Post-process results
inference.model <- postprocess.result(inference.model)
tf.inference <- inference.model$result.info$tf.inference
```

For more detailed examples, check the example scripts in the `examples` directory.
