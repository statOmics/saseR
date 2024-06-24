# saseR

saseR is a highly performant and fast algorithm to perform aberrant expression and splicing analyses. By using adapted offsets in the negative binomial framework, we have shown that aberrant and differential splicing can be modelled. Next, by using an adapted mean-variance structure, we get unbiased estimators that can be used when many latent confounders are included, which is often the case for rare Mendelian studies. With this combination of adapted offsets and fast parameter estimation, we have shown that we are much faster compared to the state-of-the-art aberrant expression and splicing algorithms, while outperforming the latter to detect aberrant splicing events.

For more information, be sure to check out our bioRxiv preprint and the vignette of this package: 
https://doi.org/10.1101/2023.06.29.547014

## NEWS

## Installation instructions

To install the development version, run;

```{r 'install_dev', eval = FALSE}
devtools::install_github("statOmics/saseR")
```

The installation should only take a few seconds.
The dependencies of the package are listed in the DESCRIPTION file of the package.

saseR is also available at Bioconductor. To install its current version, run;

```{r 'install_bioconductor', eval = FALSE}

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("saseR")

```

## Issues and bug reports

Please use https://github.com/statOmics/saseR/issues to submit issues, bug reports, and comments.

# Example data

Example data was obtained from the ASpli package.

Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz, ASpli: integrative analysis of splicing landscapes through RNA-Seq assays, Bioinformatics, Volume 37, Issue 17, September 2021, Pages 2609–2616, https://doi.org/10.1093/bioinformatics/btab141
