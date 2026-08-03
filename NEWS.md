# eDNA2IQI 4.0.0
A Non-parametric Multidimensional Scaling (NMDS) based test to check if data is
similar to training dataset used for model. If data fails test, results should
be used with caution and re-analysis or fauna analysis considered. By default,
samples failing NMDS test will return `NA` IQI results. This behavior is a
'breaking change' and version number has been bump to 4.0.0 to indicate results
may now in some cases be different to previously returned results.

# eDNA2IQI 3.0.0
eDNA2IQI now includes and adjusted IQI prediction alongside the previous prediction.  
This adjustment is correcting the predicted IQI values using the linear equation from the out-of-bag predicted IQI 
to actual IQI for the training dataset.

# eDNA2IQI 2.5.0
Random Forest Models updated to include all data from BactMetBar 1, and data provided by SEPA for the 2025 BactMetBar project extension.

# eDNA2IQI 2.4.8
Migrating to SAMS github. 
minor changes in comments and comment organization. 
locale "C" sorting added to taxa after allocation, repairing 
testthat function stability with running in an actual R session.

# eDNA2IQI 2.4.7
prediction from multiple models stable

# eDNA2IQI 2.4.6
Test Checks bugs fixed, issue with empty DF in step 3 fixed 

# eDNA2IQI 2.4.5
Test Checks added by TF, thank you!

# eDNA2IQI 2.4.4
Metadata in .csv output from prediction. 
more useful data in console at each step

# eDNA2IQI 2.4.3
added vignettes, license.md added to Rbuildignore.
example data added into inst/extdata (TEMPORARY)


# eDNA2IQI 2.4.2
external download bug fixed. note mode = "wb"
Renv snapshot Updated

# eDNA2IQI 2.4.1
Added barplot function seperately
changed examples
external download still bugged

# eDNA2IQI 2.4.0
Welcome Linux and MacOS users!
Stable for Windows and Linux, with an updated readme. 
Untested for Mac so far though. 

# eDNA2IQI 2.3.4
download check for reference taxa database
error handling if thredds server is down. tries from the Zenodo server. 

# eDNA2IQI 2.3.3
folder path doesnt need /

# eDNA2IQI 2.3.2
wrapper function for all steps

# eDNA2IQI 2.3.1
bug fix

# eDNA2IQI 2.3.0
Increase of outputs, including barplots and rarefied dataset.  

# eDNA2IQI 2.2.4
seed set at beginning and before rarefaction for iqi prediction. 
moved step 0.5 into step 1, so executes as default. 
roxygen updated

# eDNA2IQI 2.2.3
datacleaningdetails.txt now includes a line of output stating package version used

# eDNA2IQI 2.2.2
new parameteroptions file with row for cutadapt filepath

# eDNA2IQI 2.2.1
removed bug of primertrimming still left in dada2

# eDNA2IQI 2.2.0
restored cutadapt, options etc. 


# eDNA2IQI 2.1.1
fixed taxa allocation version snafu, 
removed cutadapt parameters from parameteroptions !!this may cause issues!!


# eDNA2IQI 2.1.0
removed step1.3, and trimming primers by length using the dada2:filterandTrim

# eDNA2IQI 2.0.5
biocmananager added, new parameteroptions

# eDNA2IQI 2.0.4
renaming functions to remove capitalisation
# eDNA2IQI 2.0.2
Major correction to taxa allocation format. Previous models will not work with newly processed data

# eDNA2IQI 2.0.0
package rejig

# eDNA2IQI 1.0.0

# eDNA2IQI 0.0.0.9000

* Created the package
