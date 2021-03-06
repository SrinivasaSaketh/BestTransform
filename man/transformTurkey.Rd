% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TransformationPrerequisites.R
\name{transformTurkey}
\alias{transformTurkey}
\title{Performsturkey power transformation on a vector}
\usage{
transformTurkey(x, start, end, int, plotit, verbose, quiet, statistic, returnLambda)
}
\arguments{
\item{x}{A vector of values}

\item{start}{The starting value of lambda to try}

\item{end}{The ending value of lambda to try}

\item{int}{The interval between lambda values to try}

\item{plotit}{If TRUE, produces plots of Shapiro-Wilks W or Anderson-Darling A vs. lambda, a histogram of transformed values, and a quantile-quantile plot of transformed values}

\item{verbose}{If TRUE, prints extra output for Shapiro-Wilks W or Anderson-Darling A vs. lambda}

\item{quiet}{If TRUE, doesn't print any output to the screen}

\item{statistic}{If 1, uses Shapiro-Wilks test. If 2, uses Anderson-Darling test}

\item{returnLambda}{If TRUE, returns only the lambda value, not the vector of transformed values}
}
\value{
The transformed vector of values.
}
\description{
Takes in a vector and conducts Tukey's Ladder of Powers on a vector of values to produce a more-normally distributed vector of values
}
