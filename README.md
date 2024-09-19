
<!-- README.md is generated from README.Rmd. Please edit that file -->

 <!-- badges: start -->
  [![R-CMD-check](https://github.com/Adwyness/edna2iqi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Adwyness/edna2iqi/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

# eDNA2IQI

eDNA2IQI allows users to predict IQI values from raw 16S rRNA sequence
reads. The full pipeline will manage, trim, denoise, and allocate taxa
to sequence reads. The IQI will then be predicted using a Random Forest
algorithm. Visit BactMetBar — Scottish Association for Marine Science,
Oban UK (sams.ac.uk) for more information.

## Installation

eDNA2IQI requires several packages to be installed via BiocManager,
follow the ? below to install.

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")

BiocManager::install(c("zlibbioc", "IRanges", "S4Vectors", "XVector", 
                       "GenomeInfoDb", "SparseArray", "DelayedArray", 
                       "S4Arrays", "Biobase", "Rhtslib", "BiocParallel", 
                       "Rsamtools", "Biostrings", "SummarizedExperiment", 
                       "GenomicRanges", "GenomicAlignments", "ShortRead", 
                       "dada2"))
}
```


Then you can install eDNA2IQI directly from GitHub with:

``` r
install.packages("devtools")
library(“devtools”)
devtools::install_github("SAMS-UK/eDNA2IQI")
```

Alternatively you can install manually from source:

- Download the zipped package of the full "SAMS-UK/eDNA2IQI"
  repository.
- Unzip the folder
- Open the eDNA2IQI.Rproj
- Enter the following:

``` r
install.packages("devtools")
library("devtools")
install()
```

Check this worked by running:

``` r
packageVersion("eDNA2IQI")
```

Successful installation indicated by correct version number:

``` r
[1] ‘2.4.4’
```

**Ensure to execute eDNA2IQI within a fresh R or R studio session after
manual installation**

## Documentation

The eDNA2IQI package has been split into logical steps to be called
upon. Providing the input data is in the correct format, each step
should follow on from the previous, until ultimately a predicted IQI
value is outputted.

## Getting Started

**Load data**

The first step is to ensure the input data is of the correct format. For
this example, we will be denoising a samples from a survey on
"FinFishFarmA", working from the folder path:

"C:/Users/AdW/IQI_generation/FinFishFarmA/"

Create a directory within the working directory called "input". This is
where your raw ".fastq.gz" files should be located.

These files should:

• Be generated following the SEPA protocol for the Use of DNA in
Monitoring Compliance with Seabed Mixing Zone Area Limits
(guidance-on-the-use-of-dna-in-monitoring-compliance-with-seabed-mixing-zone-area-limits.pdf
(sepa.org.uk)).

- Be paired end FASTQ files.

- Be compressed by gzip (“.gz”).

- Follow the naming standard outputted by an Illumina miseq.Eg:

>       Sample1_S31_L001_R1_001.fastq.gz
>       Sample1_S31_L001_R2_001.fastq.gz
>       Sample2_S32_L001_R1_001.fastq.gz
>       Sample2_S32_L001_R2_001.fastq.gz

Do not worry about:

- Combining samples from more than 1 sequencing run. They cannot be
  denoised together, but the package will split them and do them
  separately and recombine them.

- Some file transfer applications are not case sensitive. If case is
  lost, do not worry, the package can deal with this. But do not
  change the filename convention, or edit the files in any way.

In a fresh R or R studio session, load the eDNA2IQI library

``` r
library(eDNA2IQI)
```

Create a value for the filepath of the input directory

``` r
folder <- "C:/Users/AdW/IQI_generation/FinFishFarmA/input/"
```

## Running for the first time?

The first time the package is run after installation or update will
require the download of: the Cutadapt v4.0
[DOI:10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200) executable file during Step 1, the SILVA
138.1 reference taxa database [DOI:10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219) during Step
2. These will be firstly attempted to be downloaded from maintained links 
on the SAMS THREDDS server, then from source. Cutadapt nor the SILVA
database are our creation, we have simply mirrored them for link
stability and ease of use.
These can be downloaded before your first run through with samples, 
it is always good to be prepared:

``` r
eDNA2IQI::downloadExternal()
```

## Operating Systems

The eDNA2IQI package was built for running R on Windows OS. 
However, eDNA2IQI will yield *identical* results on Linux or Mac OS,
it just needs a little extra setup. 

Follow the installation instructions above, and also run the downloadExternal() 
above in a fresh R session.

The cutadapt.exe file is a windows executable, but can be run through R 
using "Wine", a compatability layer that needs to be installed on your OS:

On your Linux or MacOS shell
``` bash
sudo apt update
sudo apt install wine
```
Then in R
``` R
file.path(find.package("eDNA2IQI"))
```
and navigate into the /extdata/ directory containing the cutadapt_v1.exe file

then in the shell
``` bash
chmod 777 cutadapt_v1.exe
```
This will allow eDNA2IQI to execute the file from within R.

For a full walkthrough of the eDNA2IQI pipeline, check out the vignette 
```r
vignette("eDNA2IQI_Vignette.Rmd", package = "eDNA2IQI")
```

