
<!-- README.md is generated from README.Rmd. Please edit that file -->
# WARNING: No model included in the download. please save one to desktop

# eDNA2IQI v2.1

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

The eDNA2IQI packages allows users to import and clean bacterial DNA
sequence reads. Taxa names are then allocated to the sequences, using a
curated SILVA v138.1 database as a reference. The IQI can then be
predicted using a trained random forest algorithm.

## Installation

You can install eDNA2IQI directly from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("SAMS-UK/eDNA2IQI")
```

## Run full package

To run the full package:

``` r
samples <- predict_iqi_from_fastq("C:/full/file/path/to/fastq/")
```

The full package pipeline takes raw fastq.gz files, denoises the data,
assigns taxa, rareifies the data, and predicts the IQI value of the
samples using a trained random forest algorithm.

The user data is assigned taxa at the same taxa level and rarefaction
rate as the data used to train the random forest. If a sample has fewer
reads than the rarefaction rate required, the IQI can not be predicted.

## Clean data and assign taxa

To denoise data and assign taxa to samples:

``` r
samples <- readFastqSequences_assignTaxa(folder = "C:/full/file/path/to/fastq/")
```

## Predict IQI of pre-cleaned, denoised and taxa allocated data

The IQI values of data that has been denoised and taxa allocated can be
predicted with:

``` r
samples <- predict_iqi(dataframe)
```

The data must be in the form of a dataframe. Each sample is represented
by a row, with taxa as the column names. The number of reads for each
taxa in a sample is entered in the dataframe. Missing taxa must be
represented by a zero (not NA/nan/NaN).
