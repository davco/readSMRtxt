\name{getSRM}
\alias{getSRM}
\alias{getSRM.default}
\alias{print.getSRM}
\title{Convert a SRM text file into a matrix.}
\description{Convert a SRM text file into a class object. This class object include: sample name, chromatographic retention time range, number of transitions and chromatographic time-intensity matrix.}
\usage{
getSRM(x)

\method{getSRM}{default}(x)
\method{print}{getSRM}(x)
}

\arguments{
  \item{x}{ a SRM text file. }
}
\value{
 A class, basically a list including elements
 \item{sample_name}{ sample name }
 \item{timestamp}{ chromatographic retention time range }
 \item{transitions}{ number of transitions }
 \item{TI_matrix}{ chromatographic time-intensity matrix }
}

\author{David Garcia}

\examples{
path_to_spectra <- system.file("extdata", "foo_spectra.txt", package = "readSRMtxt")

spectra_class_format <- getSRM(path_to_spectra) # convert the .txt file into a class object format

# Plot a TIC
spectra_matrix <- spectra_class_format$TI_matrix
TIC <- apply(spectra_matrix[,(2:dim(spectra_matrix)[2])], 1, sum, na.rm=T)
plot(TIC~spectra_matrix[,1], type= "l", xlab="Time (min)", main=paste("TIC of", spectra_class_format$sample_name, "sample", sep=" ")) # plotting TIC
}

\keyword{convert}

