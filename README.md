
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eDNA2IQI

eDNA2IQI allows users to predict IQI values from raw 16S rRNA sequence
reads. The full pipeline will manage, trim, denoise, and allocate taxa
to sequence reads. The IQI will then be predicted using a Random Forest
algorithm. Visit BactMetBar — Scottish Association for Marine Science,
Oban UK (sams.ac.uk) for more information.

## Installation

First install several packages through BiocManager

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")

BiocManager::install(c("zlibbioc", "IRanges", "S4Vectors", "XVector", 
                       "GenomeInfoDb", "SparseArray", "DelayedArray", 
                       "S4Arrays", "Biobase", "Rhtslib", "BiocParallel", 
                       "Rsamtools", "Biostrings", "SummarizedExperiment", 
                       "GenomicRanges", "GenomicAlignments", "ShortRead", 
                       "dada2"))
```


Then you can install eDNA2IQI directly from GitHub with:

``` r
install.packages("devtools")
library(“devtools”)
devtools::install_github("SAMS-UK/eDNA2IQI")
```

Alternatively you can install manually from source:

> - Download the zipped package of the full “SAMS-UK/eDNA2IQI”
>   repository.
> - Unzip the folder
> - Open the eDNA2IQI.Rproj
> - Enter the following:

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
[1] ‘2.2.2’
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
“FinFishFarmA”, working from the folder path:

“C:/Users/AdW/IQI_generation/FinFishFarmA/”

Create a directory within the working directory called “input”. This is
where your raw “.fastq.gz” files should be located.

These files should:

• Be generated following the SEPA protocol for the Use of DNA in
Monitoring Compliance with Seabed Mixing Zone Area Limits
(guidance-on-the-use-of-dna-in-monitoring-compliance-with-seabed-mixing-zone-area-limits.pdf
(sepa.org.uk)).

> - Be paired end FASTQ files.

> - Be compressed by gzip (“.gz”).

> - Follow the naming standard outputted by an Illumina miseq.Eg:

> >       Sample1_S31_L001_R1_001.fastq.gz
> >       Sample1_S31_L001_R2_001.fastq.gz
> >       Sample2_S32_L001_R1_001.fastq.gz
> >       Sample2_S32_L001_R2_001.fastq.gz

Do not worry about:

> - Combining samples from more than 1 sequencing run. They cannot be
>   denoised together, but the package will split them and do them
>   separately and recombine them.

> - Some file transfer applications are not case sensitive. If case is
>   lost, do not worry, the package can deal with this. But do not
>   change the filename convention, or edit the files in any way.

In a fresh R or R studio session, load the eDNA2IQI library

``` r
library(eDNA2IQI)
```

Create a value for the filepath of the input directory, ///remember the forward slash at the end!///

``` r
folder <- “C:/Users/AdW/IQI_generation/FinFishFarmA/input/”
```

## Running for the first time?

The first time the package is run after installation or update will
require the download of: the cutadapt v4.0
\[<DOI:10.14806/ej.17.1.200>\] executable file during Step 1, the SILVA
138.1 reference taxa database \[<DOI:10.1093/nar/gks1219>\] during Step
2, and the most up-to-date Random Forest Model for IQI prediction. These
will be downloaded from the SAMS THREDDS server. Cutadapt nor the SILVA
database are our creation, we have simply mirrored them for link
stability and ease of use.

## Full Tutorial

**STEP 1**

This step will take your raw FASTQ files, and turn them into a
data.frame called “myASVs”, containing Amplicon Sequence Variants (ASVs)
and their abundance in each sample, through a process referred to as
denoising. It will create an output directory within the input directory
(“C:/Users/AdW/FinFishFarmA/input/outputData”), within which will be
outputted:

> - A summary of read quality for perusal and investigation following
>   any issues in IQI generation.

> - Individual batch files as .rda and .csv files.

> - dataCleaningDetails.txt, that contains all parameters used for the
>   denoising, and statistics of how many reads were lost at each stage
>   for investigation upon IQI generation issues.

> - A directory containing the intermediary files
>   (“C:/Users/AdW/FinFishFarmA/input/outputData/intermediate”). If the
>   process were to be interrupted during the denoising, the Step 1
>   function will scan these folders and resume the pipeline where it
>   left off.

Step 1 takes the folder you have previously designated, and a list of
the parameteroptions. The default for these parameters should not be
edited or amended, as this would violate the IQI prediction. However,
the quality check can be turned off by changing the value of the
“qualitycheck” parameter to “FALSE”.

``` r
data(“parameteroptions”)
parameteroptions[3,3] <- “FALSE”
```

You are ready to execute Step 1, this process can take several hours, if
in doubt that the pipeline is progressing, the intermediate directories
can be checked via windows explorer for ongoing population. Note however
that during certain steps such as error learning, no intermediate files
will be created.

``` r
myASVs <- eDNA2IQI::step1.0_readFastq(folder, parameteroptions)
```

**STEP 2**

Step 2 annotates the ASVs generated in Step 1 to “Family” level,
creating a data.frame of bacterial families, and their abundance in each
sample. The input for this step is the “myASVs” object from Step 1. ASVs
will be annotated against the SILVA 138.1 database which will be stored
locally after being downloaded on the first use.

``` r
myTAXA <- eDNA2IQI::step2.0_annotate_ASVs(folder, parameteroptions, myASVs)
```

**STEP 3**

Step 3 will take the annotated taxa data.frame from Step 2, and use the
abundances of families contained to predict an IQI for each sample. The
most up-to-date Random Forest model for IQI prediction will be
downloaded prior to this step on the first use. If a sample has fewer
reads than the number required, the IQI can not be predicted.

``` r
predicted_IQIs <- eDNA2IQI::step3.0_predict_iqi(myTAXA)
```

Predicted IQIs can then be exported as a .csv file:

``` r
write.csv(file = "predicted_IQIs.csv", predicted_IQIs)
```
