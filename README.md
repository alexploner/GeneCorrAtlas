## GeneCorrAtlas

[Bulik-Sullivan and co-authors](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3406.html) presented a matrix of genome-wide 
genetic correlations derived via [linkage-disequilibrium regression](http://www.nature.com/ng/journal/v47/n3/full/ng.3211.html). 
The resulting correlations are available both from the supplementary material of the original article and 
[from github](https://raw.githubusercontent.com/bulik/gencor/master/all_rg.csv) in an easily readable, alas somewhat dry text format 
(.csv). 

This R package wraps the data and offers simple but convenient tools for extracting  and displaying these genetic correlations 
from the R commandline.

## Installation

To install `GeneCorrAtlas` from github, you need the devtools package in your R environment. If you have not done so yet, start a fresh R session and run
```
install.packages("devtools")
```
This will install the developer tools and their dependencies. Now you can easily install `VEGAStools` via
```
devtools::install_github("alexploner/GeneCorrAtlas")
```
and you are good to go. 

## Getting started

While you are at the R command line, have a look at the package vignette:
```
vignette("GeneCorrAtlas_inR")
```
This will give you a quick overview of how to access and display the data.
