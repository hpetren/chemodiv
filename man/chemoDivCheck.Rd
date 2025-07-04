% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemoDivCheck.R
\name{chemoDivCheck}
\alias{chemoDivCheck}
\title{Check data formatting}
\usage{
chemoDivCheck(sampleData, compoundData)
}
\arguments{
\item{sampleData}{Data frame with the relative concentration of each
compound (column) in every sample (row).}

\item{compoundData}{Data frame with the compounds in \code{sampleData}
as rows. Should have a column named "compound" with common names of
the compounds, a column named "smiles" with SMILES IDs of the compounds,
and a column named "inchikey" with the InChIKey IDs for the compounds.}
}
\value{
One or several messages pointing out problems with data formatting,
or a message informing that the datasets appear to be correctly formatted.
}
\description{
Function to check that the datasets used by other functions in
the \emph{chemodiv} package are correctly formatted.
}
\details{
The function performs a number of checks on the two main datasets used
as input data, to make sure datasets are formatted in a way suitable
for the other functions in the package. This should make it easier for
users to correctly construct datasets before starting with analyses.

Two datasets are needed to use the full set of analyses included in
the package, and these can be checked for formatting issues.
The first dataset should contain data on the proportions
of different compounds (columns) in different samples (rows).
Note that all calculations of diversity, and most calculations of
dissimilarity, are only performed on relative, rather than absolute,
values. The second dataset should contain, in each of three
columns in a data frame, the compound name, SMILES and InChIKey IDs of
all the compounds present in the first dataset. See
\code{\link{chemodiv}} for details on obtaining SMILES and InChIKey IDs.
For compound names, avoid starting with a number, avoid using Greek letters
or special characters, and ensure there are no trailing spaces.
}
\examples{
data(minimalSampData)
data(minimalCompData)
chemoDivCheck(minimalSampData, minimalCompData) # Correct format
chemoDivCheck(minimalSampData, minimalCompData[c(2,3,1),]) # Incorrect format

data(alpinaSampData)
data(alpinaCompData)
chemoDivCheck(sampleData = alpinaSampData, compoundData = alpinaCompData)
}
