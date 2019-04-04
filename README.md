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
library("readSRMtxt")
path_to_spectra <- system.file("extdata", "foo_spectra.txt", package = "readSRMtxt")

     spectra_class_format <- getSRM(path_to_spectra) # convert the .txt file into a class object format

     # Plot a TIC
     spectra_matrix <- spectra_class_format$TI_matrix
     TIC <- apply(spectra_matrix[,(2:dim(spectra_matrix)[2])], 1, sum, na.rm=T)
     plot(TIC~spectra_matrix[,1], type= "l", xlab="Time (min)", main=paste("TIC of", spectra_class_format$sample_name, "sample", sep=" ")) # plotting TIC
```
