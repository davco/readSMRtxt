# readSRMtxt


# Description
Convert MRM and Schedule-MRM text file into a class object. This class object include: sample name, chromatographic retention time range, number of transitions and chromatographic time-intensity matrix (TI_matrix). TI_matrix contain: a time variable in the first column, intensity signal variables in the following columns for each transition.

# Install
```r
library(devtools)
install_github("davco/readSRMtxt")
```

# Install

1. Convert .wiff into a text file using ProteoWizard Software (http://proteowizard.sourceforge.net/) 
2. Convert text file into a time-intensity matrix for each transition with this package.


```r
library(devtools)
install_github("davco/readSRMtxt")
How to convert .wiff file into a matrix
```
