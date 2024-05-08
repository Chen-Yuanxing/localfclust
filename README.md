# localfclust

The package **localfclust** implements the local functional clustering (LFC) method proposed by Chen et al. (2024). LFC aims to identify both the global clustering structure and multiple different local clustering structures.
The method relies on the locally group-wise fusion and B-spline approximation.

In addition to the developed package, we also provide some files to reproduce the simulation results and analysis results of the COVID-19 data.
Specifically, the file *simulation_replicate.R* implements the LFC method based on the settings of Example 1, and evaluates global/local clustering accuracy and estimation accuracy for the obtained results.
Another file *covid19_analysis.R* applies the LFC method to the COVID-19 data, stored in the Rdata object *covid19_normalized_nytimes.RData*,
and outputs multiple figures, which correspond to Figures 1, 5, and 6 of Chen et al. (2024).

## Installation

The development version can be installed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Chen-Yuanxing/localfclust")
```
