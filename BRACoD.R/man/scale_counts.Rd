% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BRACoD.R
\name{scale_counts}
\alias{scale_counts}
\title{Normalize OTU counts and add a pseudo count}
\usage{
scale_counts(df_counts)
}
\arguments{
\item{df_counts}{A dataframe of OTU counts. Samples are rows and bacteria are columns.}
}
\value{
a dataframe of relative abundance data
}
\description{
BRACoD requires relative abundance and cannot handle zeros, so this function
adds a small pseudo count (1/10th the smallest non-zero value).
}
