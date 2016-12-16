# *flippant*

An `R` package for the automated analysis of dithionite scramblase assays.


## Installation

To install the development version, you first need the *devtools* package.

```{r}
install.packages("devtools")
```

Then you can install the *flippant* package using

```{r}
library(devtools)
install_bitbucket("jog2030/flippant")
```

A stable version is planned to be made available via Bioconductor.

# The manuscript

The manuscript is provided in `Rmd` form in the `inst/manuscript` directory.  Note that building the manuscript takes around an hour, since it runs the whole analysis.

# The dataset

The original data in contained in a `zip` file in the `inst/extdata` directory.  Each file within the `zip` archive is a text file.  Lines 4 through to the penultimate line contain tab delimited text data.

# Functions

`scramblase_assay_traces` draws time series of relative fluorescence vs. acquisition time, possibly split by experiment.

`scramblase_assay_plot` draws plots of the probability of at least one flippase vs. the protein to phospholipid ratio.

`scramblase_assay_stats` calculates the fit constant for each experiment.

`scramblase_assay_calculations` provides a lower-level overview of model fit details.