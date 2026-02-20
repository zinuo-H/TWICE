
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TWICE

<!-- badges: start -->

<!-- badges: end -->

A New Enrichment Method for Exposomics: Topology-Weighted Integrated
Correlation Enrichment (TWICE) Analysis

The goal of TWiCE is to combine the key strengths of existing
approaches: 1. Effectively utilize correlation values, similar to
Functional Class Scoring (FCS) methods, so that genes with stronger
correlations contribute more to the pathway score. 2. Incorporate
pathway topology information, as in Pathway Topology-based (PT) methods,
to assign greater weights to biologically important genes such as
upstream or highly connected nodes. In this way, TWiCE bridges the gap
between correlation analysis and functional enrichment, providing a more
biologically meaningful interpretation of exposome-omics data.

## Installation

You can install the released version of `TWICE` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("TWICE")
```

You can install the development version of `TWICE` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Asa12138/pcutils")
devtools::install_github("Asa12138/TWICE")
```

## Example

This is a basic example:

``` r
library(TWICE)
## basic example code
load_all_pathway_NS() -> all_pathway_NS
#> ==============================load all_pathway_NS===============================
#> ================all_pathway_NS time: 2026-02-20 21:52:33.167294=================
#> If you want to update all_pathway_NS, use `update_all_pathway_NS()`
# Generate an example data box
cor_df <- data.frame(
  name = all_pathway_NS$name,
  cor = runif(nrow(all_pathway_NS), 0, 1)
)
cor_df <- dplyr::sample_n(cor_df, 200)

calculate_TWICE(all_pathway_NS, cor_df, mode = ) -> res

head(res)
#>   Pathway_id K_num Exist_K_num
#> 1    ko05216     2           2
#> 2    ko05332     3           3
#> 3    ko03263     1           1
#> 4    ko05033     1           1
#> 5    ko04740     2           2
#> 6    ko04930     5           5
#>                                             Exist_K     S_mean      S_sum
#> 1                               ko:K09289|ko:K04371 0.04322915 0.08645830
#> 2                     ko:K03156|ko:K05405|ko:K04519 0.03456180 0.10368539
#> 3                                         ko:K25295 0.06200105 0.06200105
#> 4                                        cpd:C01996 0.06213583 0.06213583
#> 5                              cpd:C00575|ko:K04515 0.03205475 0.06410950
#> 6 ko:K04527|ko:K03156|ko:K04440|ko:K07203|ko:K04371 0.02205663 0.11028314
#>   p_value p_adjust Direction
#> 1   0.001   0.1805  Positive
#> 2   0.001   0.1805  Positive
#> 3   0.002   0.1805  Positive
#> 4   0.002   0.1805  Positive
#> 5   0.003   0.1805  Positive
#> 6   0.003   0.1805  Positive
```
