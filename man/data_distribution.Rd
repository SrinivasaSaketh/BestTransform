% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataDistribution.R
\name{data_distribution}
\alias{data_distribution}
\title{Find the distribution of a dataframe}
\usage{
data_distribution(data, "dv")
}
\arguments{
\item{data}{A data frame for which the distribution needs to be identified}

\item{dv}{A string value specifying the column name of DV}
}
\value{
A data frame with different distribution metrics such as variable type, is_normal, shapiro p value, skewness score, kurtosis score, pearson score, etc.
}
\description{
Takes in a data frame and identifies the data distribution
}
